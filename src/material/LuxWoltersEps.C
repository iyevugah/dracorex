//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LuxWoltersEps.h"

registerMooseObject("dracorexApp", LuxWoltersEps);
registerMooseObject("dracorexApp", ADLuxWoltersEps);

template <bool is_ad>
InputParameters
LuxWoltersEpsTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnCreepStressUpdateBaseTempl<is_ad>::validParams();
  params.addClassDescription(
      "This class uses the stress update material in a radial return isotropic creep model"
      "This class computes the modified Lubby2 creep.");
  // Maxwell parameters
  params.addParam<Real>("mvM", "Maxwell viscosity parameter");
  params.addRequiredParam<Real>("etaM0", "Initial Maxwell Viscosity");
  // Kelvin parameters
  params.addParam<Real>("mvK", "Kelvin ViscoParameter");
  params.addParam<Real>("mk", "Kelvin Elastic Parameter");
  params.addRequiredParam<Real>("etaK0", "Initial Kelvin Viscosity");
  params.addParam<Real>("GK0", "Initial Kelvin Shear Modulus");
  //Model Params
  params.addRequiredParam<Real>("a4", "model parameter 4");
  params.addRequiredParam<Real>("a5", "model parameter 5");
  params.addRequiredParam<Real>("a6", "model parameter 6");
  params.addRequiredParam<Real>("a7", "model parameter 7");
  params.addRequiredParam<Real>("a8", "model parameter 8");
  params.addRequiredParam<Real>("a9", "model parameter 9");
  params.addRequiredParam<Real>("a10", "model parameter 10");
  params.addRequiredParam<Real>("a15", "model parameter 15");
  params.addRequiredParam<Real>("a16", "model parameter 16");
  params.addRequiredParam<Real>("a17", "model parameter 17");
  params.addRequiredParam<Real>("sigma0", "model parameter sigma0");
  params.addRequiredParam<Real>("L", "model parameter L");
  params.addRequiredParam<Real>("L1", "model parameter L1");
  params.addRequiredParam<Real>("T", "model parameter T");
  params.addRequiredParam<Real>("a", "model parameter a");
  params.addRequiredParam<Real>("b", "model parameter b");
  return params;
}

template <bool is_ad>
LuxWoltersEpsTempl<is_ad>::LuxWoltersEpsTempl(const InputParameters & parameters)
  : RadialReturnCreepStressUpdateBaseTempl<is_ad>(parameters),
    _etaM0(this->template getParam<Real>("etaM0")),
    _mvM(this->template getParam<Real>("mvM")),
    _mvK(this->template getParam<Real>("mvK")),
    _mk(this->template getParam<Real>("mk")),
    _etaK0(this->template getParam<Real>("etaK0")),
    _GK0(this->template getParam<Real>("GK0")),
    _kelvin_creep_rate(this->template declareGenericProperty<Real, is_ad>("kelvin_creep_rate")),
    _kelvin_creep_rate_old(this->template getMaterialPropertyOld<Real>("kelvin_creep_rate")),
    // damage parameter
    _damage_param(this->template declareGenericProperty<Real, is_ad>(this->_base_name +"damage_param")),
    _damage_param_old(this->template getMaterialPropertyOld<Real>(this->_base_name + "damage_param")),
    // Model Params
    _a4 (this->template getParam<Real>("a4")),
    _a5 (this->template getParam<Real>("a5")),
    _a6 (this->template getParam<Real>("a6")),
    _a7 (this->template getParam<Real>("a7")),
    _a8 (this->template getParam<Real>("a8")),
    _a9 (this->template getParam<Real>("a9")),
    _a10 (this->template getParam<Real>("a10")),
    _a15 (this->template getParam<Real>("a15")),
    _a16 (this->template getParam<Real>("a16")),
    _a17 (this->template getParam<Real>("a17")),
    _sigma0 (this->template getParam<Real>("sigma0")),
    _L (this->template getParam<Real>("L")),
    _L1 (this->template getParam<Real>("L1")),
    _T (this->template getParam<Real>("T")),
    _a (this->template getParam<Real>("a")),
    _b (this->template getParam<Real>("b"))
{
  if (_etaM0 == 0.0 && _etaK0 == 0.0)
    mooseError("LuxWoltersEps: at least one of the creep should be active.");
}

template <bool is_ad>
void
LuxWoltersEpsTempl<is_ad>::initQpStatefulProperties()
{
  _kelvin_creep_rate[_qp] = 0.0;
  _damage_param[_qp] = 0.0; // Initialize damage parameter
}

template <bool is_ad>
void
LuxWoltersEpsTempl<is_ad>::propagateQpStatefulProperties()
{
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp];
  _damage_param[_qp] = _damage_param_old[_qp];
  RadialReturnStressUpdateTempl<is_ad>::propagateQpStatefulPropertiesRadialReturn();
}


template <bool is_ad>
void
LuxWoltersEpsTempl<is_ad>::computeStressInitialize(const GenericReal<is_ad> & /*effective_trial_stress*/,
                                               const GenericRankFourTensor<is_ad> & /*elasticity_tensor*/,
                                          RankTwoTensor & stress_new, const RankTwoTensor deviatoric_stress)
{
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp];
  _damage_param[_qp] = _damage_param_old[_qp];
  Real J3_sig = deviatoric_stress.det(); // 3rd invariant of deviatoric stress
  Real smin = ((stress_new) - deviatoric_stress).trace(); // minimum principal stress
}

/// Compute the Residuals when automatic_differentiation = false
template <bool is_ad>
template <typename ScalarType>
ScalarType
LuxWoltersEpsTempl<is_ad>::computeResidualInternal(const GenericReal<is_ad> & effective_trial_stress,
                                                 const ScalarType & scalar)
{
  const ScalarType stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  //computing the Damage Parameter
  ScalarType lode_param = (27 / 2) * (J3_sig) / pow(stress_delta,3.0);  // Lode Parameter
  ScalarType lode_ang = 1 - ((2 / pi) * acos(lode_param));       // Lode Angle
  ScalarType eta_D = (1.0) - (_a4 * std::exp(-_a5*smin));
  ScalarType beta_TC = (_a6) - (_a7* std::exp(-_a8*smin));
  ScalarType kappa_beta = pow(1.0 / (cos(lode_ang + (pi / 6.0)) + _a9 * sin(lode_ang + (pi / 6.0))),
                                exp(-(_a10) * smin));
  ScalarType F_ds = stress_delta - ((eta_D)*(beta_TC)*(kappa_beta));
  ScalarType F_dz = 6.0 * (-smin);
  _damage_param[_qp] = _damage_param_old[_qp];
  const ScalarType  damage_rate = (_a15/(pow(1.0 - _damage_param[_qp],_a17 ))) * pow((F_ds + F_dz),_a16);
  _damage_param[_qp] = _damage_param_old[_qp] + MetaPhysicL::raw_value(damage_rate) * _dt; // update damage param

  const ScalarType etaM = _etaM0 * pow((stress_delta/_sigma0)/(1-_damage_param[_qp]),_a)
                          * std::exp(_mvM * stress_delta)*std::exp(_L*_T);
  const ScalarType etaK = _etaK0 * std::exp(_mvK * stress_delta);
  const ScalarType GK = _GK0 * pow((stress_delta/_sigma0)/(1-_damage_param[_qp]),_b)
                          * std::exp(_mk * stress_delta)*std::exp(_L1*_T);
  _kelvin_creep_rate[_qp] =  _kelvin_creep_rate_old[_qp] + MetaPhysicL::raw_value(scalar);

  if (_etaM0 != 0.0 && _etaK0 != 0.0)
      {
      // Maxwell and Kelvin
       const ScalarType M_creep_rate = (stress_delta / ((1.0 -_damage_param[_qp]) * etaM));
       const ScalarType K_creep_rate = (stress_delta/((1.0-_damage_param[_qp])* etaK)) -
                                   ((GK*(_kelvin_creep_rate[_qp]))/ ((std::sqrt(2.0/3.0))* etaK));
       return (M_creep_rate + K_creep_rate) * _dt - scalar;
      }

       else if (_etaM0 != 0.0 && _etaK0 == 0.0)
      // Maxwell
      {
        const ScalarType creep_rate = (stress_delta / ((1.0 -_damage_param[_qp]) * etaM));
        return creep_rate * _dt - scalar;
      }
      // Kelvin
       const ScalarType creep_rate = (stress_delta/((1.0-_damage_param[_qp])* etaK)) -
                                ((GK*(_kelvin_creep_rate[_qp]))/ ((std::sqrt(2.0/3.0))* etaK));
         return creep_rate * _dt - scalar;
 }


 /// Compute the Derivatives when automatic_differentiation = false
 /// in the SingleVariableReturnMappingSolution.C file in MOOSE. Otherwise,
 /// no need to compute (or use this) derivatives.
 template <bool is_ad>
 GenericReal<is_ad>
 LuxWoltersEpsTempl<is_ad>::computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                             const GenericReal<is_ad> & scalar)
 {
   const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
   // computing the Damage Parameter
    GenericReal<is_ad> lode_param = (27 / 2) * (J3_sig) / (stress_delta);  // Lode Parameter
   const GenericReal<is_ad> lode_ang = 1 - ((2 / pi) * acos(lode_param));       // Lode Angle
   const GenericReal<is_ad> eta_D = (1.0) - (_a4 * std::exp(-_a5*smin));
   const GenericReal<is_ad> beta_TC = (_a6) - (_a7* std::exp(-_a8*smin));
   const GenericReal<is_ad> kappa_beta = pow(1.0 / (cos(lode_ang + (pi / 6.0)) + _a9 * sin(lode_ang + (pi / 6.0))),
                                 exp(-(_a10) * smin));
   const GenericReal<is_ad> F_ds = stress_delta - ((eta_D)*(beta_TC)*(kappa_beta));
   const GenericReal<is_ad> F_dz = 6.0 * (-smin);
   _damage_param[_qp] = _damage_param_old[_qp];
   const GenericReal<is_ad> damage_rate = (_a15/(1.0 - pow(_damage_param[_qp],_a17 ))) * pow((F_ds + F_dz),_a16);
   _damage_param[_qp] = _damage_param_old[_qp] + damage_rate * _dt;    // update damage param

   const GenericReal<is_ad> etaM = _etaM0 * pow((stress_delta/_sigma0)/(1-_damage_param[_qp]),_a)
                                  * std::exp(_mvM * stress_delta)*std::exp(_L*_T);
   const GenericReal<is_ad> etaK = _etaK0 * std::exp(_mvK * stress_delta);
   const GenericReal<is_ad> GK = _GK0 * pow((stress_delta/_sigma0)/(1-_damage_param[_qp]),_b)
                                  * std::exp(_mk * stress_delta)*std::exp(_L1*_T);
   _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp] + scalar;

   if (_etaM0 != 0.0 && _etaK0 != 0.0)
   {
     const GenericReal<is_ad> M_creep_rate_derivative =
       (_three_shear_modulus/((1.0 -_damage_param[_qp]) * etaM)) * (-1+(stress_delta*_mvM));
   // Use the following for the M_creep_rate derivative if GK must be differentiated:
   // const GenericReal<is_ad> M_creep_rate_derivative = (_three_shear_modulus/((1.0 -_damage_param[_qp]) * etaM)) *
   //                                                (-1 + ((stress_delta*_mvM ) - (_a*_damage_param[_qp])+_a));
     const GenericReal<is_ad> K_creep_rate_derivative =
       (_three_shear_modulus/etaK)* (((stress_delta*_mvK)/(1.0 -_damage_param[_qp])) - (1.0/(1.0-_damage_param[_qp])) -
       ((_kelvin_creep_rate[_qp]*GK*_mvK)/(std::sqrt(2.0/3.0))));

   // Use the following for the K_creep_rate derivative if GK must be differentiated:
   // const GenericReal<is_ad> K_creep_rate_derivative =
   //  (_three_shear_modulus/etaK)* (((stress_delta*_mvK)/(1.0 -_damage_param[_qp])) - (1.0/(1.0-_damage_param[_qp])) -
   //   ((_kelvin_creep_rate[_qp]*GK*_mvK)/(std::sqrt(2.0/3.0))) +
   //     (_kelvin_creep_rate[_qp]*GK*(_mvK*stress_delta - (_b*_damage_param[_qp]) +(_b))));

 return (M_creep_rate_derivative + K_creep_rate_derivative) * _dt - 1.0;
  }

    else if (_etaM0 != 0.0 && _etaK0 == 0.0)
{
  // Maxwell
  const GenericReal<is_ad> creep_rate_derivative =
        (_three_shear_modulus/((1.0 -_damage_param[_qp]) * etaM)) * ((stress_delta*_mvM)-1);
  // Use the following for the M_creep_rate derivative if GK must be differentiated:
  // const GenericReal<is_ad> M_creep_rate_derivative = (_three_shear_modulus/((1.0 -_damage_param[_qp]) * etaM)) *
  //                                                (-1 + ((stress_delta*_mvM ) - (_a*_damage_param[_qp])+_a)); 
  return creep_rate_derivative * _dt - 1.0;
 }
 //Kelvin
 const GenericReal<is_ad> creep_rate_derivative =
 (_three_shear_modulus/etaK)* (((stress_delta*_mvK)/(1.0 -_damage_param[_qp])) - (1.0/(1.0-_damage_param[_qp])) -
 ((_kelvin_creep_rate[_qp]*GK*_mvK)/(std::sqrt(2.0/3.0))));
 // Use the following for the Kelvin creep_rate derivative if GK must be differentiated:
 // const GenericReal<is_ad> creep_rate_derivative =
 //  (_three_shear_modulus/etaK)* (((stress_delta*_mvK)/(1.0 -_damage_param[_qp])) - (1.0/(1.0-_damage_param[_qp])) -
 //   ((_kelvin_creep_rate[_qp]*GK*_mvK)/(std::sqrt(2.0/3.0))) +
 //     (_kelvin_creep_rate[_qp]*GK*(_mvK*stress_delta - (_b*_damage_param[_qp]) +(_b))));
return creep_rate_derivative * _dt - 1.0;
}


template <bool is_ad>
void
LuxWoltersEpsTempl<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _creep_strain[_qp] += plastic_strain_increment;
}

template <bool is_ad>
void
LuxWoltersEpsTempl<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <bool is_ad>
bool
LuxWoltersEpsTempl<is_ad>::substeppingCapabilityEnabled()
{
  return this->template getParam<bool>("use_substep");
}

template class LuxWoltersEpsTempl<false>;
template class LuxWoltersEpsTempl<true>;

template Real LuxWoltersEpsTempl<false>::computeResidualInternal<Real>(const Real &, const Real &);
template ADReal LuxWoltersEpsTempl<true>::computeResidualInternal<ADReal>(const ADReal &,
                                                                     const ADReal &);
template ChainedReal
LuxWoltersEpsTempl<false>::computeResidualInternal<ChainedReal>(const Real &, const ChainedReal &);

template ChainedADReal
LuxWoltersEpsTempl<true>::computeResidualInternal<ChainedADReal>(const ADReal &, const ChainedADReal &);
