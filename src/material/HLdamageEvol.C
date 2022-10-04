#include "HLdamageEvol.h"

registerMooseObject("TensorMechanicsApp", HLdamageEvol);
registerMooseObject("TensorMechanicsApp", ADHLdamageEvol);

template <bool is_ad>
InputParameters
HLdamageEvolTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnCreepStressUpdateBaseTempl<is_ad>::validParams();
  params.addClassDescription("This class takes advantage of the newton iteration"
     "in RadialReturnStressUpdate to solve for the incremental value of a damage "
     "variable. Based on this incremental value, the evolution of the damage"
     "variable is then solved iteratively.");

  // model and input parameters
  params.addRequiredParam<Real>("a4", "model parameter 4");
  params.addRequiredParam<Real>("a5", "model parameter 5");
  params.addRequiredParam<Real>("a6", "model parameter 6");
  params.addRequiredParam<Real>("a7", "model parameter 7");
  params.addRequiredParam<Real>("a8", "model parameter 8");
  params.addRequiredParam<Real>("a9", "model parameter 9");
  params.addRequiredParam<Real>("a10", "model paramete 10");
  params.addRequiredParam<Real>("a15", "model paramete 15");
  params.addRequiredParam<Real>("a17", "model paramete 17");
//  params.addParam<MaterialPropertyName>("damage_property", "HLdamage",
//                                      "The name of the damage property.");
//  params.addRequiredParam<std::string>("effective_inelastic_strain_name",
//  "Name of the material property that stores the effective inelastic strain");
  return params;
}


template <bool is_ad>
HLdamageEvolTempl<is_ad>::HLdamageEvolTempl(const InputParameters & parameters)
  : RadialReturnCreepStressUpdateBaseTempl<is_ad>(parameters),
    _a4 (this->template getParam<Real>("a4")),
    _a5 (this->template getParam<Real>("a5")),
    _a6 (this->template getParam<Real>("a6")),
    _a7 (this->template getParam<Real>("a7")),
    _a8 (this->template getParam<Real>("a8")),
    _a9 (this->template getParam<Real>("a9")),
    _a10 (this->template getParam<Real>("a10")),
    _a15 (this->template getParam<Real>("a15")),
   _a17 (this->template getParam<Real>("a17")),
//damage_variable(this->template declareGenericProperty<Real, is_ad>(
//this->_base_name +
//this->template getParam<std::string>("effective_inelastic_strain_name"))),
//damage_variable_old(this->template getMaterialPropertyOld<Real>(
//this->_base_name +
//this->template getParam<std::string>("effective_inelastic_strain_name")))
damage_variable(this->template declareGenericProperty<Real, is_ad>("damage_index_val")),
damage_variable_old(this->template getMaterialPropertyOld<Real>("damage_index_val"))

{
}


template <bool is_ad>
void
HLdamageEvolTempl<is_ad>::initQpStatefulProperties()
{
  damage_variable[_qp] = 0.0;
}


template <bool is_ad>
void
HLdamageEvolTempl<is_ad>::propagateQpStatefulProperties()
{
  damage_variable[_qp] = damage_variable_old[_qp];
}


 template <bool is_ad>
 void
 HLdamageEvolTempl<is_ad>::computeStressInitialize (const GenericReal<is_ad> & /*effective_trial_stress*/,
                                                    const GenericRankFourTensor<is_ad> & /*elasticity_tensor*/,
                                                    RankTwoTensor & stress_new, const RankTwoTensor deviatoric_stress)
{
  Real J3_sig = deviatoric_stress.det(); // 3rd invariant of deviatoric stress
  Real smin = ((stress_new) - deviatoric_stress).trace(); // minimum principal stress
 }

template <bool is_ad>
template <typename ScalarType>
ScalarType
HLdamageEvolTempl<is_ad>::computeResidualInternal(const GenericReal<is_ad> & effective_trial_stress,
                                              const ScalarType & scalar)
{
  const ScalarType stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const ScalarType lode_param = (27 / 2) * (J3_sig) / (stress_delta);  // Lode Parameter
  const ScalarType lode_ang = 1 - ((2 / pi) * acos(lode_param));       // Lode Angle
  const ScalarType eta_D = (1.0) - (_a4 * std::exp(-_a5*smin));
  const ScalarType beta_TC = (_a6) - (_a7* std::exp(-_a8*smin));
  const ScalarType kappa_beta = pow(1.0 / (cos(lode_ang + (pi / 6.0)) + _a9 * sin(lode_ang + (pi / 6.0))),
                                exp(-(_a10) * smin));
  const ScalarType F_ds = stress_delta - ((eta_D)*(beta_TC)*(kappa_beta));
  const ScalarType F_dz = 6.0 * (-smin);
  const ScalarType damage_rate = (_a15/(1.0-(damage_rate * _dt))) * (F_ds + F_dz);
  //const ScalarType _scalar = (damage_rate * _dt);
  return (damage_rate * _dt) - scalar;

}

template <bool is_ad>
GenericReal<is_ad>
HLdamageEvolTempl<is_ad>::computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                        const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
//  const GenericReal<is_ad> _scalar = scalar;
  const GenericReal<is_ad> damage_rate_derivative =
      -(_three_shear_modulus);
  return damage_rate_derivative * _dt - 1.0;
}



template <bool is_ad>
void
HLdamageEvolTempl<is_ad>::iterationFinalize(GenericReal<is_ad> scalar)
{
  damage_variable[_qp] =
    damage_variable_old[_qp] + scalar;
}


/*
template <bool is_ad>
void
HLdamageEvolTempl<is_ad>::computeQpProperties()
{
damage_variable[_qp] =
  damage_variable_old[_qp] + _scalar;
}
*/


/*
template <bool is_ad>
void
HLdamageEvolTempl<is_ad>::computedamagestate(const GenericReal<is_ad> & scalar)
{

  damage_variable[_qp] =
    damage_variable_old[_qp] + scalar;
}
*/


/*
template <bool is_ad>
void
HLdamageEvolTempl<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <bool is_ad>
bool
HLdamageEvolTempl<is_ad>::substeppingCapabilityEnabled()
{
  return this->template getParam<bool>("use_substep");
}
*/


template class HLdamageEvolTempl<false>;
template class HLdamageEvolTempl<true>;

template Real HLdamageEvolTempl<false>::computeResidualInternal<Real>(const Real &, const Real &);
template ADReal HLdamageEvolTempl<true>::computeResidualInternal<ADReal>(const ADReal &,
                                                                     const ADReal &);
template ChainedReal
HLdamageEvolTempl<false>::computeResidualInternal<ChainedReal>(const Real &, const ChainedReal &);

template ChainedADReal
HLdamageEvolTempl<true>::computeResidualInternal<ChainedADReal>(const ADReal &, const ChainedADReal &);
