//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EpsMmodLubby2.h"

registerMooseObject("TensorMechanicsApp", EpsMmodLubby2);
registerMooseObject("TensorMechanicsApp", ADEpsMmodLubby2);

template <bool is_ad>
InputParameters
EpsMmodLubby2Templ<is_ad>::validParams()
{
  InputParameters params = RadialReturnCreepStressUpdateBaseTempl<is_ad>::validParams();
  params.addClassDescription("This class uses the stress update material in a radial return isotropic creep model"
    "This class computes the secondary creep component of the modified Lubby2 creep.");


  // Linear strain hardening parameters
  params.addRequiredParam<Real>("mvM", "Maxwell viscosity parameter");
  params.addRequiredParam<Real>("etaM0", "Initial Maxwell Viscosity");
  return params;
}

template <bool is_ad>
EpsMmodLubby2Templ<is_ad>::EpsMmodLubby2Templ(const InputParameters & parameters)
  : RadialReturnCreepStressUpdateBaseTempl<is_ad>(parameters),
    _mvM(this->template getParam<Real>("mvM")),
    _etaM0(this->template getParam<Real>("etaM0"))
{
}


template <bool is_ad>
template <typename ScalarType>
ScalarType
EpsMmodLubby2Templ<is_ad>::computeResidualInternal(const GenericReal<is_ad> & effective_trial_stress,
                                              const ScalarType & scalar)
{
  const ScalarType stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const ScalarType etaMax = _etaM0 * std::exp(_mvM * stress_delta);
  const ScalarType creep_rate = stress_delta / (3.0 * etaMax);

  return creep_rate * _dt - scalar;
}


template <bool is_ad>
GenericReal<is_ad>
EpsMmodLubby2Templ<is_ad>::computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                        const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const GenericReal<is_ad> etaMax = _etaM0 * std::exp(_mvM * effective_trial_stress);
  const GenericReal<is_ad> creep_rate_derivative = -_three_shear_modulus / (3.0 * etaMax) * (1.0 + _mvM * stress_delta);

  return creep_rate_derivative * _dt - 1.0;
}


template <bool is_ad>
void
EpsMmodLubby2Templ<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _creep_strain[_qp] += plastic_strain_increment;
}

template <bool is_ad>
void
EpsMmodLubby2Templ<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <bool is_ad>
bool
EpsMmodLubby2Templ<is_ad>::substeppingCapabilityEnabled()
{
  return this->template getParam<bool>("use_substep");
}

template class EpsMmodLubby2Templ<false>;
template class EpsMmodLubby2Templ<true>;

template Real EpsMmodLubby2Templ<false>::computeResidualInternal<Real>(const Real &, const Real &);
template ADReal EpsMmodLubby2Templ<true>::computeResidualInternal<ADReal>(const ADReal &,
                                                                     const ADReal &);
template ChainedReal
EpsMmodLubby2Templ<false>::computeResidualInternal<ChainedReal>(const Real &, const ChainedReal &);

template ChainedADReal
EpsMmodLubby2Templ<true>::computeResidualInternal<ChainedADReal>(const ADReal &, const ChainedADReal &);
