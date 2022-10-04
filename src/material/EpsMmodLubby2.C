//* AJ

#include "EpsMmodLubby2.h"

#include "Function.h"
#include "ElasticityTensorTools.h"

registerMooseObject("TensorMechanicsApp", ADEpsMmodLubby2);
registerMooseObject("TensorMechanicsApp", EpsMmodLubby2);

template <bool is_ad>
InputParameters
EpsMmodLubby2Templ<is_ad>::validParams()
{
  InputParameters params = RadialReturnStressUpdateTempl<is_ad>::validParams();
  params.addClassDescription("This class uses the stress update material in a radial return isotropic creep model"
    "This class computes the primary creep component of the modified Lubby2 creep.");

  // input parameters
  params.addRequiredParam<Real>("mvM", "Maxwell viscosity parameter");
  params.addRequiredParam<Real>("etaM0", "Initial Maxwell Viscosity");
  return params;
}

template <bool is_ad>
EpsMmodLubby2Templ<is_ad>::EpsMmodLubby2Templ(
    const InputParameters & parameters)
  : RadialReturnStressUpdateTempl<is_ad>(parameters),
    _plastic_strain(this->template declareGenericProperty<RankTwoTensor, is_ad>("plastic_strain")),
    _plastic_strain_old(this->template getMaterialPropertyOld<RankTwoTensor>("plastic_strain")),
    _mvM(this->template getParam<Real>("mvM")),
    _etaM0(this->template getParam<Real>("etaM0"))
{
}


template <bool is_ad>
void
EpsMmodLubby2Templ<is_ad>::initQpStatefulProperties()
{
  _plastic_strain[_qp].zero();
}

template <bool is_ad>
void
EpsMmodLubby2Templ<is_ad>::propagateQpStatefulProperties()
{
  _plastic_strain[_qp] = _plastic_strain_old[_qp];

  EpsMmodLubby2Templ<is_ad>::propagateQpStatefulPropertiesRadialReturn();
}

template <bool is_ad>
void
EpsMmodLubby2Templ<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & effective_trial_stress,
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  _plastic_strain[_qp] = _plastic_strain_old[_qp];
}

template <bool is_ad>
GenericReal<is_ad>
EpsMmodLubby2Templ<is_ad>::computeResidual(
    const GenericReal<is_ad> & effective_trial_stress, const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const GenericReal<is_ad> etaMax = _etaM0 * std::exp(_mvM * stress_delta);
  const GenericReal<is_ad> creep_rate = stress_delta / (3.0 * etaMax);
  return creep_rate * _dt - scalar;
}

template <bool is_ad>
GenericReal<is_ad>
EpsMmodLubby2Templ<is_ad>::computeDerivative(
    const GenericReal<is_ad> & effective_trial_stress, const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const GenericReal<is_ad> etaMax = _etaM0 * std::exp(_mvM * effective_trial_stress);
  const GenericReal<is_ad> creep_rate_derivative = _three_shear_modulus / (3.0 * etaMax) * (1.0 + _mvM * stress_delta);
  return creep_rate_derivative * _dt - 1.0;
}


template <bool is_ad>
void
EpsMmodLubby2Templ<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _plastic_strain[_qp] += plastic_strain_increment;
}


template class EpsMmodLubby2Templ<false>;
template class EpsMmodLubby2Templ<true>;
