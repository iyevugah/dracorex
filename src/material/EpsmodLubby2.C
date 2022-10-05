//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EpsmodLubby2.h"

registerMooseObject("TensorMechanicsApp", EpsmodLubby2);
registerMooseObject("TensorMechanicsApp", ADEpsmodLubby2);

template <bool is_ad>
InputParameters
EpsmodLubby2Templ<is_ad>::validParams()
{
  InputParameters params = RadialReturnCreepStressUpdateBaseTempl<is_ad>::validParams();
  params.addClassDescription(
      "This class uses the stress update material in a radial return isotropic creep model"
      "This class computes the modified Lubby2 creep.");
  // Maxwell parameters
  params.addParam<Real>("mvM", "Maxwell viscosity parameter");
  params.addParam<Real>("n_exponent", "Exponent on effective stress in power-law equation");
  params.addRequiredParam<Real>("etaM0", "Initial Maxwell Viscosity");
  // Kelvin parameters
  params.addParam<Real>("mvK", "Kelvin ViscoParameter");
  params.addParam<Real>("mk", "Kelvin Elastic Parameter");
  params.addRequiredParam<Real>("etaK0", "Initial Kelvin Viscosity");
  params.addParam<Real>("GK0", "Initial Kelvin Shear Modulus");
  return params;
}

template <bool is_ad>
EpsmodLubby2Templ<is_ad>::EpsmodLubby2Templ(const InputParameters & parameters)
  : RadialReturnCreepStressUpdateBaseTempl<is_ad>(parameters),
    _etaM0(this->template getParam<Real>("etaM0")),
    _mvM(this->template getParam<Real>("mvM")),
    _n_exponent(this->template getParam<Real>("n_exponent")),
    _mvK(this->template getParam<Real>("mvK")),
    _mk(this->template getParam<Real>("mk")),
    _etaK0(this->template getParam<Real>("etaK0")),
    _GK0(this->template getParam<Real>("GK0")),
    _kelvin_creep_rate(this->template declareGenericProperty<Real, is_ad>("kelvin_creep_rate")),
    _kelvin_creep_rate_old(this->template getMaterialPropertyOld<Real>("kelvin_creep_rate"))
{
  if (_etaM0 == 0.0 && _etaK0 == 0.0)
    mooseError("EpsmodLubby2: at least one of the creep should be active.");
}

template <bool is_ad>
void
EpsmodLubby2Templ<is_ad>::initQpStatefulProperties()
{
  _kelvin_creep_rate[_qp] = 0.0;
}

template <bool is_ad>
void
EpsmodLubby2Templ<is_ad>::propagateQpStatefulProperties()
{
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp];
  RadialReturnStressUpdateTempl<is_ad>::propagateQpStatefulPropertiesRadialReturn();
}

template <bool is_ad>
void
EpsmodLubby2Templ<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & effective_trial_stress,
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp];
}

template <bool is_ad>
GenericReal<is_ad>
EpsmodLubby2Templ<is_ad>::computeResidual(const GenericReal<is_ad> & effective_trial_stress,
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
EpsmodLubby2Templ<is_ad>::computeResidualMK(const GenericReal<is_ad> & effective_trial_stress,
                                            const GenericReal<is_ad> & scalar)
{
  GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  // Maxwell
  const GenericReal<is_ad> etaM = _etaM0 * std::exp(_mvM * stress_delta);
  const GenericReal<is_ad> M_creep_rate = stress_delta / (3.0 * etaM);
  // Kelvin
  const GenericReal<is_ad> etaK = _etaK0 * std::exp(_mvK * stress_delta);
  const GenericReal<is_ad> GK = _GK0 * std::exp(_mk * stress_delta);
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp] + scalar;
  const GenericReal<is_ad> creep_rate = stress_delta / (3.0 * etaK) - (_three_shear_modulus*_kelvin_creep_rate[_qp])/(3.0*etaK);
  return creep_rate * _dt - scalar;
}

template <bool is_ad>
GenericReal<is_ad>
EpsmodLubby2Templ<is_ad>::computeResidualM(const GenericReal<is_ad> & effective_trial_stress,
                                           const GenericReal<is_ad> & scalar)
{
  GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  // Maxwell
  const GenericReal<is_ad> etaM = _etaM0 * std::exp(_mvM * stress_delta);
  const GenericReal<is_ad> creep_rate = stress_delta / (3.0 * etaM);
  return creep_rate * _dt - scalar;
}

template <bool is_ad>
GenericReal<is_ad>
EpsmodLubby2Templ<is_ad>::computeResidualK(const GenericReal<is_ad> & effective_trial_stress,
                                           const GenericReal<is_ad> & scalar)
{
  GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const GenericReal<is_ad> etaK = _etaK0 * std::exp(_mvK * stress_delta);
  const GenericReal<is_ad> GK = _GK0 * std::exp(_mk * stress_delta);
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp] + scalar;
  const GenericReal<is_ad> creep_rate = stress_delta / (3.0 * etaK) - (_three_shear_modulus*_kelvin_creep_rate[_qp])/(3.0*etaK);
  return creep_rate * _dt - scalar;
}

template <bool is_ad>
GenericReal<is_ad>
EpsmodLubby2Templ<is_ad>::computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
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
EpsmodLubby2Templ<is_ad>::computeDerivativeMK(const GenericReal<is_ad> & effective_trial_stress,
                                              const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  // Maxwell
  const GenericReal<is_ad> etaM = _etaM0 * std::exp(_mvM * effective_trial_stress);
  const GenericReal<is_ad> M_creep_rate_derivative =
      _three_shear_modulus / (3.0 * etaM) * (1.0 + _mvM * stress_delta);
  // Kelvin
  const GenericReal<is_ad> etaK = _etaK0 * std::exp(_mvK * stress_delta);
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp] + scalar;
  const GenericReal<is_ad> K_creep_rate_derivative =
      (_three_shear_modulus / (3.0 * etaK)) * (1.0 + (_mvK * (stress_delta - (_three_shear_modulus*_kelvin_creep_rate[_qp]))));

  return (M_creep_rate_derivative + K_creep_rate_derivative) * _dt - 1.0;
}

template <bool is_ad>
GenericReal<is_ad>
EpsmodLubby2Templ<is_ad>::computeDerivativeM(const GenericReal<is_ad> & effective_trial_stress,
                                             const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const GenericReal<is_ad> etaM = _etaM0 * std::exp(_mvM * stress_delta);
  const GenericReal<is_ad> creep_rate_derivative =
      _three_shear_modulus / (3.0 * etaM) * (1.0 + _mvM * stress_delta);

  return creep_rate_derivative * _dt - 1.0;
}

template <bool is_ad>
GenericReal<is_ad>
EpsmodLubby2Templ<is_ad>::computeDerivativeK(const GenericReal<is_ad> & effective_trial_stress,
                                             const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  // Kelvin
  const GenericReal<is_ad> etaK = _etaK0 * std::exp(_mvK * stress_delta);
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp] + scalar;
  const GenericReal<is_ad> creep_rate_derivative =
      (_three_shear_modulus / (3.0 * etaK)) * (1.0 + (_mvK * (stress_delta - (_three_shear_modulus*_kelvin_creep_rate[_qp]))));

  return creep_rate_derivative * _dt - 1.0;
}

template <bool is_ad>
Real
EpsmodLubby2Templ<is_ad>::computeStrainEnergyRateDensity(
    const GenericMaterialProperty<RankTwoTensor, is_ad> & stress,
    const GenericMaterialProperty<RankTwoTensor, is_ad> & strain_rate)
{
  if (_n_exponent <= 1)
    return 0.0;

  Real creep_factor = _n_exponent / (_n_exponent + 1);

  return MetaPhysicL::raw_value(creep_factor * stress[_qp].doubleContraction((strain_rate)[_qp]));
}

template <bool is_ad>
void
EpsmodLubby2Templ<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _creep_strain[_qp] += plastic_strain_increment;
}

template <bool is_ad>
void
EpsmodLubby2Templ<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <bool is_ad>
bool
EpsmodLubby2Templ<is_ad>::substeppingCapabilityEnabled()
{
  return this->template getParam<bool>("use_substep");
}

template class EpsmodLubby2Templ<false>;
template class EpsmodLubby2Templ<true>;

// template Real EpsmodLubby2Templ<false>::computeResidualInternal<Real>(const Real &, const Real
// &); template ADReal EpsmodLubby2Templ<true>::computeResidualInternal<ADReal>(const ADReal &,
//                                                                          const ADReal &);
// template ChainedReal
// EpsmodLubby2Templ<false>::computeResidualInternal<ChainedReal>(const Real &, const ChainedReal
// &);
//
// template ChainedADReal
// EpsmodLubby2Templ<true>::computeResidualInternal<ChainedADReal>(const ADReal &,
//                                                                 const ChainedADReal &);
