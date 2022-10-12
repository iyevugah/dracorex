//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EpsKmodLubby2.h"

registerMooseObject("TensorMechanicsApp", EpsKmodLubby2);
registerMooseObject("TensorMechanicsApp", ADEpsKmodLubby2);

template <bool is_ad>
InputParameters
EpsKmodLubby2Templ<is_ad>::validParams()
{
  InputParameters params = RadialReturnCreepStressUpdateBaseTempl<is_ad>::validParams();
  params.addClassDescription("This class uses the stress update material in a radial return isotropic creep model"
    "This class computes the secondary creep component of the modified Lubby2 creep.");


  // Linear strain hardening parameters
  params.addRequiredParam<Real>("mvK", "Kelvin ViscoParameter");
  params.addRequiredParam<Real>("mk", "Kelvin Elastic Parameter");
  params.addRequiredParam<Real>("etaK0", "Initial Kelvin Viscosity");
  params.addRequiredParam<Real>("GK0", "Initial Kelvin Shear Modulus");
  return params;
}

template <bool is_ad>
EpsKmodLubby2Templ<is_ad>::EpsKmodLubby2Templ(const InputParameters & parameters)
  : RadialReturnCreepStressUpdateBaseTempl<is_ad>(parameters),
  _effective_creep_increment(this->template declareGenericProperty<Real, is_ad>("effective_creep_increment")),
  _effective_creep_increment_old(this->template getMaterialPropertyOld<Real>("effective_creep_increment")),
  _mvK(this->template getParam<Real>("mvK")),
  _mk(this->template getParam<Real>("mk")),
  _etaK0 (this->template getParam<Real>("etaK0")),
  _GK0(this->template getParam<Real>("GK0"))
{
}

template <bool is_ad>
void
EpsKmodLubby2Templ<is_ad>::initQpStatefulProperties()
{
  _effective_creep_increment[_qp] = 0.0;
}

template <bool is_ad>
void
EpsKmodLubby2Templ<is_ad>::propagateQpStatefulProperties()
{
  _effective_creep_increment[_qp] = _effective_creep_increment_old[_qp];

  RadialReturnStressUpdateTempl<is_ad>::propagateQpStatefulPropertiesRadialReturn();
}

template <bool is_ad>
template <typename ScalarType>
ScalarType
EpsKmodLubby2Templ<is_ad>::computeResidualInternal(const GenericReal<is_ad> & effective_trial_stress,
                                              const ScalarType & scalar)
{
  const ScalarType stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const ScalarType etaKMax = _etaK0 * std::exp(_mvK * stress_delta);
  const ScalarType GKMax = _GK0 * std::exp(_mk * stress_delta);
  _effective_creep_increment[_qp] = _effective_creep_increment_old[_qp] + MetaPhysicL::raw_value(scalar);
  const ScalarType creep_rate = stress_delta / (3.0 * etaKMax) - (GKMax*_effective_creep_increment[_qp])/(etaKMax);

  return creep_rate * _dt - scalar;
}


template <bool is_ad>
GenericReal<is_ad>
EpsKmodLubby2Templ<is_ad>::computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                        const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const GenericReal<is_ad> etaKMax = _etaK0 * std::exp(_mvK * stress_delta);
  const GenericReal<is_ad> GKMax = _GK0 * std::exp(_mk * stress_delta);
  _effective_creep_increment[_qp] = _effective_creep_increment_old[_qp] + scalar;
  const GenericReal<is_ad> creep_rate_derivative = (_three_shear_modulus / (3.0 * etaKMax)) * (1.0 + (_mvK * (stress_delta - (_three_shear_modulus*_effective_creep_increment[_qp]))));

  return creep_rate_derivative * _dt - 1.0;
}


template <bool is_ad>
void
EpsKmodLubby2Templ<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _creep_strain[_qp] += plastic_strain_increment;
}

template <bool is_ad>
void
EpsKmodLubby2Templ<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <bool is_ad>
bool
EpsKmodLubby2Templ<is_ad>::substeppingCapabilityEnabled()
{
  return this->template getParam<bool>("use_substep");
}

template class EpsKmodLubby2Templ<false>;
template class EpsKmodLubby2Templ<true>;

template Real EpsKmodLubby2Templ<false>::computeResidualInternal<Real>(const Real &, const Real &);
template ADReal EpsKmodLubby2Templ<true>::computeResidualInternal<ADReal>(const ADReal &,
                                                                     const ADReal &);
template ChainedReal
EpsKmodLubby2Templ<false>::computeResidualInternal<ChainedReal>(const Real &, const ChainedReal &);

template ChainedADReal
EpsKmodLubby2Templ<true>::computeResidualInternal<ChainedADReal>(const ADReal &, const ChainedADReal &);
