#include "Lubby2.h"
#include <RankTwoTensor.h>

registerMooseObject("dracorexApp", Lubby2);
registerMooseObject("dracorexApp", ADLubby2);

template <bool is_ad>
InputParameters
Lubby2Templ<is_ad>::validParams()
{
  InputParameters params = RadialReturnCreepStressUpdateBaseTempl<is_ad>::validParams();
  params.addClassDescription(
    "This is a basic creep model for rocksalt based on the Lubby2 model. See Zhang"
    "and Nagel (2020): Error-controlled implicit time integration of elasto-visco-plastic"
    "constitutive models for rock salt. It uses the StressUpdate class in the radialreturn"
    "isotropic creep model");
  // Maxwell parameters
  params.addParam<Real>("mvM", "Maxwell viscosity parameter");
  params.addRequiredParam<Real>("etaM0", "Initial Maxwell Viscosity");
  // Kelvin parameters
  params.addParam<Real>("mvK", "Kelvin ViscoParameter");
  params.addParam<Real>("mk", "Kelvin Elastic Parameter");
  params.addRequiredParam<Real>("etaK0", "Initial Kelvin Viscosity");
  params.addParam<Real>("GK0", "Initial Kelvin Shear Modulus");
  params.addParam<RealTensorValue>("EQstress", "Applied Equivalent Stress");
//  params.addParam<FunctionName>("EQstressFn","The applied equivalent stress as a function");
//  params.addParam<bool>("function_EQstress",
//                    false,
//                    "Whether the applied stress is a single value or a function");
 return params;
}

template <bool is_ad>
Lubby2Templ<is_ad>::Lubby2Templ(const InputParameters & parameters)
  : RadialReturnCreepStressUpdateBaseTempl<is_ad>(parameters),
    _etaM0(this->template getParam<Real>("etaM0")),
    _mvM(this->template getParam<Real>("mvM")),
    _mvK(this->template getParam<Real>("mvK")),
    _mk(this->template getParam<Real>("mk")),
    _etaK0(this->template getParam<Real>("etaK0")),
    _GK0(this->template getParam<Real>("GK0")),
    _kelvin_creep_rate(this->template declareGenericProperty<Real, is_ad>("kelvin_creep_rate")),
    _kelvin_creep_rate_old(this->template getMaterialPropertyOld<Real>("kelvin_creep_rate")),
//    _EQstress(this->template getParam<Real>("EQstress")),
//    _EQstressFn(isParamValid("EQstressFn") ? &getFunction("_EQstressFn"): nullptr),
//    _function_EQstress(parameters.get<bool>("function_EQstress"))
    _EQstress(parameters.isParamValid("EQstress")
              ? this->template getParam<RealTensorValue>("EQstress")
              : RealTensorValue(5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)) //RankTwoTensor
{
  if (_etaM0 == 0.0 && _etaK0 == 0.0)
    mooseError("Lubby2: at least one of the creep should be active.");
}



template <bool is_ad>
void
Lubby2Templ<is_ad>::initQpStatefulProperties()
{
  _kelvin_creep_rate[_qp] = 0.0;
}



template <bool is_ad>
void
Lubby2Templ<is_ad>::propagateQpStatefulProperties()
{
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp];
  RadialReturnStressUpdateTempl<is_ad>::propagateQpStatefulPropertiesRadialReturn();
}



template <bool is_ad>
void
Lubby2Templ<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & effective_trial_stress,
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp];
}



/// Compute the Residuals when automatic_differentiation = false
template <bool is_ad>
template <typename ScalarType>
ScalarType
Lubby2Templ<is_ad>::computeResidualInternal(const GenericReal<is_ad> & effective_trial_stress,
                                                 const ScalarType & scalar)
{
   RankTwoTensor _sigma;
//   if (_function_EQstress)
   // Use the applied stress function from the input file.
//   {
//      _sigma =_function_EQstress;
//   }
//   else
     // Use the (single value) applied stress from the input file.
// {
      _sigma = _EQstress;
// }

  RankTwoTensor deviator_stress = _sigma.deviatoric();
  // compute the J2 stress
  ScalarType dev_stress_squared =
      deviator_stress.doubleContraction(deviator_stress);
  ScalarType J2 = dev_stress_squared/2.0;
  ScalarType eff = std::sqrt(3*J2);

  const ScalarType stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);

  const ScalarType etaM = _etaM0 * std::exp(_mvM * eff);
  const ScalarType etaK = _etaK0 * std::exp(_mvK * eff);
  const ScalarType GK = _GK0 * std::exp(_mk * eff);
  _kelvin_creep_rate[_qp] =  _kelvin_creep_rate_old[_qp] + MetaPhysicL::raw_value(scalar);
//  _kelvin_creep_rate[_qp] += MetaPhysicL::raw_value(scalar);


    if (_etaM0 != 0.0 && _etaK0 != 0.0)
  {
    // Maxwell and Kelvin
  const ScalarType M_creep_rate = stress_delta / (3.0 * etaM);

//  ScalarType K_creep_rate; //const ScalarType
//  _kelvin_creep_rate[_qp] += MetaPhysicL::raw_value(K_creep_rate * _dt);
 const ScalarType K_creep_rate = (stress_delta / (3.0 * etaK)) - ((GK *_kelvin_creep_rate[_qp])/(etaK));
//  const ScalarType K_creep_rate = (stress_delta / (3.0 * etaK)) - ((GK*_kelvin_creep_rate[_qp]*std::sqrt(2./3.))/(etaK));

return ((M_creep_rate + K_creep_rate) * _dt) - scalar;

   }
    else if (_etaM0 != 0.0 && _etaK0 == 0.0)
  // Maxwell
  {
   const ScalarType creep_rate = stress_delta / (3.0 * etaM);

return (creep_rate * _dt) - scalar;

  }
   else
  // Kelvin
  {
//  ScalarType creep_rate; //const ScalarType
//  _kelvin_creep_rate[_qp] += MetaPhysicL::raw_value(creep_rate * _dt);

  const ScalarType creep_rate =
      //(stress_delta / (3.0 * etaK)) - ((GK*_kelvin_creep_rate[_qp]*std::sqrt(2./3.))/etaK);
 (stress_delta / (3.0 * etaK)) - ((GK*_kelvin_creep_rate[_qp])/(etaK));

 return (creep_rate * _dt) - scalar;
  }
}


/// Compute the Derivatives when automatic_differentiation = false
template <bool is_ad>
GenericReal<is_ad>
Lubby2Templ<is_ad>::computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                            const GenericReal<is_ad> & scalar)
{
/*
  RankTwoTensor _sigma;
//  if (_function_EQstress)
  // Use the applied stress function from the input file.
//  {
//     _sigma =_function_EQstress;
//  }
//    else
    // Use the (single value) applied stress from the input file.
//      {
        _sigma = _EQstress;
//      }

  RankTwoTensor deviator_stress = _sigma.deviatoric();
   // compute the J2 stress
  GenericReal<is_ad> dev_stress_squared =
    deviator_stress.doubleContraction(deviator_stress);
  GenericReal<is_ad> J2 = dev_stress_squared/2.0;
  GenericReal<is_ad> eff = std::sqrt(3*J2);

  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const GenericReal<is_ad> etaM = _etaM0 * std::exp(_mvM * eff);
  const GenericReal<is_ad> etaK = _etaK0 * std::exp(_mvK * eff);
  const GenericReal<is_ad> GK = _GK0 * std::exp(_mk * eff);
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp] + scalar;

  if (_etaM0 != 0.0 && _etaK0 != 0.0)
   // Maxwell and Kelvin
  {
    const GenericReal<is_ad> M_creep_rate_derivative =
      -_three_shear_modulus / (3.0 * etaM) * (1.0 + (_mvM * stress_delta));
  const GenericReal<is_ad> K_creep_rate_derivative =
      (_three_shear_modulus / (3.0 * etaK)) * (-1.0 + (_mvK * (stress_delta - (3.0 *GK*_kelvin_creep_rate[_qp]))) + (3.0 *GK*_kelvin_creep_rate[_qp]*_mk));

return (M_creep_rate_derivative + K_creep_rate_derivative) * _dt - 1.0;
  }
      else if (_etaM0 != 0.0 && _etaK0 == 0.0)
    // Maxwell
  {
    const GenericReal<is_ad> creep_rate_derivative =
          -_three_shear_modulus / (3.0 * etaM) * (1.0 + (_mvM * stress_delta));
    return creep_rate_derivative * _dt - 1.0;
   }
   else
   // Kelvin
   {
  const GenericReal<is_ad> creep_rate_derivative =
        (_three_shear_modulus / (3.0 * etaK)) * (-1.0 + (_mvK * (stress_delta - (3.0 *GK*_kelvin_creep_rate[_qp]))) + (3.0 *GK*_kelvin_creep_rate[_qp]*_mk));

 return creep_rate_derivative * _dt - 1.0;
   }
*/
  return 0.0;
}



template <bool is_ad>
void
Lubby2Templ<is_ad>::computeStressFinalize(
const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _creep_strain[_qp] += plastic_strain_increment;
}

template <bool is_ad>
void
Lubby2Templ<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <bool is_ad>
bool
Lubby2Templ<is_ad>::substeppingCapabilityEnabled()
{
  return this->template getParam<bool>("use_substep");
}

template class Lubby2Templ<false>;
template class Lubby2Templ<true>;

template Real Lubby2Templ<false>::computeResidualInternal<Real>(const Real &, const Real &);
template ADReal Lubby2Templ<true>::computeResidualInternal<ADReal>(const ADReal &,
                                                                     const ADReal &);
template ChainedReal
Lubby2Templ<false>::computeResidualInternal<ChainedReal>(const Real &, const ChainedReal &);

template ChainedADReal
Lubby2Templ<true>::computeResidualInternal<ChainedADReal>(const ADReal &, const ChainedADReal &);
