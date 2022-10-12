//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EpsKHouLux.h"

registerMooseObject("TensorMechanicsApp", EpsKHouLux);
registerMooseObject("TensorMechanicsApp", ADEpsKHouLux);

template <bool is_ad>
InputParameters
EpsKHouLuxTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnCreepStressUpdateBaseTempl<is_ad>::validParams();
  params.addClassDescription("This class uses the stress update material in a radial return isotropic creep model"
    "This class computes the secondary creep component of the modified Lubby2 creep.");


  /// Kelvin parameters
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
  params.addRequiredParam<Real>("a10", "model paramete 10");
  params.addRequiredParam<Real>("a15", "model paramete 15");
  params.addRequiredParam<Real>("a17", "model paramete 17");
  return params;
}

template <bool is_ad>
EpsKHouLuxTempl<is_ad>::EpsKHouLuxTempl(const InputParameters & parameters)
  : RadialReturnCreepStressUpdateBaseTempl<is_ad>(parameters),
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
   _a17 (this->template getParam<Real>("a17"))
{
}


template <bool is_ad>
void
EpsKHouLuxTempl<is_ad>::initQpStatefulProperties()
{
  _kelvin_creep_rate[_qp] = 0.0;
  _damage_param[_qp] = 0.0; // Initialize damage parameter
}

template <bool is_ad>
void
EpsKHouLuxTempl<is_ad>::propagateQpStatefulProperties()
{
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp];
  _damage_param[_qp] = _damage_param_old[_qp];
  RadialReturnStressUpdateTempl<is_ad>::propagateQpStatefulPropertiesRadialReturn();
}

template <bool is_ad>
void
EpsKHouLuxTempl<is_ad>::computeStressInitialize(const GenericReal<is_ad> & /*effective_trial_stress*/,
                                               const GenericRankFourTensor<is_ad> & /*elasticity_tensor*/,
                                          RankTwoTensor & stress_new, const RankTwoTensor deviatoric_stress)
{
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp];
  _damage_param[_qp] = _damage_param_old[_qp];
  Real J3_sig = deviatoric_stress.det(); // 3rd invariant of deviatoric stress
  Real smin = ((stress_new) - deviatoric_stress).trace(); // minimum principal stress
}

/// Damage evolution
template <bool is_ad>
void
EpsKHouLuxTempl<is_ad>::updateDamageParam(const GenericReal<is_ad> & effective_trial_stress,
                                          const GenericReal<is_ad> & scalar)
{
  GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  //Damage Parameter
  const GenericReal<is_ad> lode_param = (27 / 2) * (J3_sig) / (stress_delta);  // Lode Parameter
  const GenericReal<is_ad> lode_ang = 1 - ((2 / pi) * acos(lode_param));       // Lode Angle
  const GenericReal<is_ad> eta_D = (1.0) - (_a4 * std::exp(-_a5*smin));
  const GenericReal<is_ad> beta_TC = (_a6) - (_a7* std::exp(-_a8*smin));
  const GenericReal<is_ad> kappa_beta = pow(1.0 / (cos(lode_ang + (pi / 6.0)) + _a9 * sin(lode_ang + (pi / 6.0))),
                                exp(-(_a10) * smin));
  const GenericReal<is_ad> F_ds = stress_delta - ((eta_D)*(beta_TC)*(kappa_beta));
  const GenericReal<is_ad> F_dz = 6.0 * (-smin);
  _damage_param[_qp] = _damage_param_old[_qp];
  const GenericReal<is_ad> damage_rate = (_a15/(1.0-(_damage_param[_qp]))) * (F_ds + F_dz);
  _damage_param[_qp] = _damage_param_old[_qp] + damage_rate * _dt;
}


/// Compute the Residuals when automatic_differentiation = false
template <bool is_ad>
template <typename ScalarType>
ScalarType
EpsKHouLuxTempl<is_ad>::computeResidualInternal(const GenericReal<is_ad> & effective_trial_stress,
                                              const ScalarType & scalar)
{
  const ScalarType stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const ScalarType etaK = _etaK0 * std::exp(_mvK * stress_delta);
  const ScalarType GK = _GK0 * std::exp(_mk * stress_delta);
  _kelvin_creep_rate[_qp] =  _kelvin_creep_rate_old[_qp] + MetaPhysicL::raw_value(scalar);
  const ScalarType creep_rate = (stress_delta/((1.0-_damage_param[_qp])* etaK)) -
                                ((GK*(_kelvin_creep_rate[_qp]))/ ((std::sqrt(2.0/3.0))* etaK));

  return creep_rate * _dt - scalar;
}


template <bool is_ad>
GenericReal<is_ad>
EpsKHouLuxTempl<is_ad>::computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                        const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const GenericReal<is_ad> etaK = _etaK0 * std::exp(_mvK * stress_delta);
  const GenericReal<is_ad> GK = _GK0 * std::exp(_mk * stress_delta);
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp] + scalar;
  const GenericReal<is_ad> creep_rate_derivative =
  (_three_shear_modulus/etaK)* (((stress_delta*_mvK)/(1.0 -_damage_param[_qp])) - (1.0/(1.0-_damage_param[_qp])) -
  ((_kelvin_creep_rate[_qp]*GK*_mvK)/(std::sqrt(2.0/3.0)))); //+ ((_kelvin_creep_rate[_qp]*GK*_mk)/(std::sqrt(2.0/3.0))) Use this term if it is OK to include derivative of GK

 return creep_rate_derivative * _dt - 1.0;
}


template <bool is_ad>
void
EpsKHouLuxTempl<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _creep_strain[_qp] += plastic_strain_increment;
}

template <bool is_ad>
void
EpsKHouLuxTempl<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <bool is_ad>
bool
EpsKHouLuxTempl<is_ad>::substeppingCapabilityEnabled()
{
  return this->template getParam<bool>("use_substep");
}

template class EpsKHouLuxTempl<false>;
template class EpsKHouLuxTempl<true>;

template Real EpsKHouLuxTempl<false>::computeResidualInternal<Real>(const Real &, const Real &);
template ADReal EpsKHouLuxTempl<true>::computeResidualInternal<ADReal>(const ADReal &,
                                                                     const ADReal &);
template ChainedReal
EpsKHouLuxTempl<false>::computeResidualInternal<ChainedReal>(const Real &, const ChainedReal &);

template ChainedADReal
EpsKHouLuxTempl<true>::computeResidualInternal<ChainedADReal>(const ADReal &, const ChainedADReal &);
