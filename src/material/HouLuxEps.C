//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HouLuxEps.h"

registerMooseObject("dracorexApp", HouLuxEps);
registerMooseObject("dracorexApp", ADHouLuxEps);

template <bool is_ad>
InputParameters
HouLuxEpsTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnCreepStressUpdateBaseTempl<is_ad>::validParams();
  params.addClassDescription(
      "This class uses the stress update material in a radial return isotropic creep model"
      "This class computes the modified Lubby2 creep.");
  // Maxwell parameters
  params.addRequiredParam<MaterialPropertyName>("damage_index",
                                          "Name of the material property containing the "
                                          "damage index, which goes from 0 (undamaged) to 1 "
                                          "(fully damaged)");
  params.addParam<Real>("mvM", "Maxwell viscosity parameter");
  params.addRequiredParam<Real>("etaM0", "Initial Maxwell Viscosity");
  // Kelvin parameters
  params.addParam<Real>("mvK", "Kelvin ViscoParameter");
  params.addParam<Real>("mk", "Kelvin Elastic Parameter");
  params.addRequiredParam<Real>("etaK0", "Initial Kelvin Viscosity");
  params.addParam<Real>("GK0", "Initial Kelvin Shear Modulus");
  return params;
}

template <bool is_ad>
HouLuxEpsTempl<is_ad>::HouLuxEpsTempl(const InputParameters & parameters)
  : RadialReturnCreepStressUpdateBaseTempl<is_ad>(parameters),
    _etaM0(this->template getParam<Real>("etaM0")),
    _mvM(this->template getParam<Real>("mvM")),
    _n_exponent(this->template getParam<Real>("n_exponent")),
    _mvK(this->template getParam<Real>("mvK")),
    _mk(this->template getParam<Real>("mk")),
    _etaK0(this->template getParam<Real>("etaK0")),
    _GK0(this->template getParam<Real>("GK0")),
    _kelvin_creep_rate(this->template declareGenericProperty<Real, is_ad>("kelvin_creep_rate")),
    _kelvin_creep_rate_old(this->template getMaterialPropertyOld<Real>("kelvin_creep_rate")),
    _damage_property(this->template getGenericMaterialProperty<Real, is_ad>("damage_index"))
{
  if (_etaM0 == 0.0 && _etaK0 == 0.0)
    mooseError("HouLuxEps: at least one of the creep should be active.");
}

template <bool is_ad>
void
HouLuxEpsTempl<is_ad>::initQpStatefulProperties()
{
  _kelvin_creep_rate[_qp] = 0.0;
}

template <bool is_ad>
void
HouLuxEpsTempl<is_ad>::propagateQpStatefulProperties()
{
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp];
  RadialReturnStressUpdateTempl<is_ad>::propagateQpStatefulPropertiesRadialReturn();
}

template <bool is_ad>
void
HouLuxEpsTempl<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & effective_trial_stress,
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp];
}

template <bool is_ad>
GenericReal<is_ad>
HouLuxEpsTempl<is_ad>::computeResidual(const GenericReal<is_ad> & effective_trial_stress,
                                          const GenericReal<is_ad> & scalar)
{
  if (_etaM0 != 0.0 && _etaK0 != 0.0)
    return computeResidualMK(effective_trial_stress, scalar);
  else if (_etaM0 != 0.0 && _etaK0 == 0.0)
    return computeResidualM(effective_trial_stress, scalar);
  return computeResidualK(effective_trial_stress, scalar);
}

template <bool is_ad>
GenericReal<is_ad>
HouLuxEpsTempl<is_ad>::computeResidualMK(const GenericReal<is_ad> & effective_trial_stress,
                                            const GenericReal<is_ad> & scalar)
{
  GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  // Maxwell
  const GenericReal<is_ad> etaM = _etaM0 * std::exp(_mvM * stress_delta);
  const GenericReal<is_ad> M_creep_rate = (stress_delta / ((1.0 -_damage_property[_qp]) * etaM));
  // Kelvin
  const GenericReal<is_ad> etaK = _etaK0 * std::exp(_mvK * stress_delta);
  const GenericReal<is_ad> GK = _GK0 * std::exp(_mk * stress_delta);
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp] + scalar;
  const GenericReal<is_ad> creep_rate = (stress_delta/((1.0-_damage_property[_qp])* etaK)) -
  ((GK*(_kelvin_creep_rate[_qp]))/ ((std::sqrt(2.0/3.0))* etaK));
  return creep_rate * _dt - scalar;
}

template <bool is_ad>
GenericReal<is_ad>
HouLuxEpsTempl<is_ad>::computeResidualM(const GenericReal<is_ad> & effective_trial_stress,
                                           const GenericReal<is_ad> & scalar)
{
  GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  // Maxwell
  const GenericReal<is_ad> etaM = _etaM0 * std::exp(_mvM * stress_delta);
  const GenericReal<is_ad> creep_rate = (stress_delta / ((1.0 -_damage_property[_qp]) * etaM));
  return creep_rate * _dt - scalar;
}

template <bool is_ad>
GenericReal<is_ad>
HouLuxEpsTempl<is_ad>::computeResidualK(const GenericReal<is_ad> & effective_trial_stress,
                                           const GenericReal<is_ad> & scalar)
{
  GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const GenericReal<is_ad> etaK = _etaK0 * std::exp(_mvK * stress_delta);
  const GenericReal<is_ad> GK = _GK0 * std::exp(_mk * stress_delta);
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp] + scalar;
  const GenericReal<is_ad> creep_rate = (stress_delta/((1.0-_damage_property[_qp])* etaK)) -
  ((GK*(_kelvin_creep_rate[_qp]))/ ((std::sqrt(2.0/3.0))* etaK));
  return creep_rate * _dt - scalar;
}

template <bool is_ad>
GenericReal<is_ad>
HouLuxEpsTempl<is_ad>::computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                            const GenericReal<is_ad> & scalar)
{
  if (_etaM0 != 0.0 && _etaK0 != 0.0)
    return computeDerivativeMK(effective_trial_stress, scalar);
  else if (_etaM0 != 0.0 && _etaK0 == 0.0)
    return computeDerivativeM(effective_trial_stress, scalar);
  return computeDerivativeK(effective_trial_stress, scalar);
}

template <bool is_ad>
GenericReal<is_ad>
HouLuxEpsTempl<is_ad>::computeDerivativeMK(const GenericReal<is_ad> & effective_trial_stress,
                                              const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  // Maxwell
  const GenericReal<is_ad> etaM = _etaM0 * std::exp(_mvM * effective_trial_stress);
  const GenericReal<is_ad> M_creep_rate_derivative =
      (_three_shear_modulus/((1.0-_damage_property[_qp]) * etaM)) * ((stress_delta*_mvM)-1);
  // Kelvin
  const GenericReal<is_ad> etaK = _etaK0 * std::exp(_mvK * stress_delta);
  const GenericReal<is_ad> GK = _GK0 * std::exp(_mk * stress_delta);
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp] + scalar;
  const GenericReal<is_ad> K_creep_rate_derivative =
  (_three_shear_modulus/etaK)* (((stress_delta*_mvK)/(1.0-_damage_property[_qp])) - (1.0/(1.0-_damage_property[_qp])) -
  ((_kelvin_creep_rate[_qp]*GK*_mvK)/(std::sqrt(2.0/3.0))) - ((_kelvin_creep_rate[_qp]*GK*_mk)/(std::sqrt(2.0/3.0))));

  return (M_creep_rate_derivative + K_creep_rate_derivative) * _dt - 1.0;
}

template <bool is_ad>
GenericReal<is_ad>
HouLuxEpsTempl<is_ad>::computeDerivativeM(const GenericReal<is_ad> & effective_trial_stress,
                                             const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const GenericReal<is_ad> etaM = _etaM0 * std::exp(_mvM * stress_delta);
  const GenericReal<is_ad> creep_rate_derivative =
      (_three_shear_modulus/((1.0-_damage_property[_qp]) * etaM)) * ((stress_delta*_mvM)-1);

  return creep_rate_derivative * _dt - 1.0;
}

template <bool is_ad>
GenericReal<is_ad>
HouLuxEpsTempl<is_ad>::computeDerivativeK(const GenericReal<is_ad> & effective_trial_stress,
                                             const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  // Kelvin
  const GenericReal<is_ad> etaK = _etaK0 * std::exp(_mvK * stress_delta);
  const GenericReal<is_ad> GK = _GK0 * std::exp(_mk * stress_delta);
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp] + scalar;
  const GenericReal<is_ad> creep_rate_derivative = (_three_shear_modulus/etaK)* (((stress_delta*_mvK)/(1.0-_damage_property[_qp])) - (1.0/(1.0-_damage_property[_qp])) -
  ((_kelvin_creep_rate[_qp]*GK*_mvK)/(std::sqrt(2.0/3.0))) - ((_kelvin_creep_rate[_qp]*GK*_mk)/(std::sqrt(2.0/3.0))));

  return creep_rate_derivative * _dt - 1.0;
}

template <bool is_ad>
void
HouLuxEpsTempl<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _creep_strain[_qp] += plastic_strain_increment;
}

template <bool is_ad>
void
HouLuxEpsTempl<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <bool is_ad>
bool
HouLuxEpsTempl<is_ad>::substeppingCapabilityEnabled()
{
  return this->template getParam<bool>("use_substep");
}

template class HouLuxEpsTempl<false>;
template class HouLuxEpsTempl<true>;
