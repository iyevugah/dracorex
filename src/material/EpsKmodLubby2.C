#include "EpsKmodLubby2.h"

#include "Function.h"
#include "ElasticityTensorTools.h"

registerMooseObject("dracorexApp", ADEpsKmodLubby2);
registerMooseObject("dracorexApp", EpsKmodLubby2);

template <bool is_ad>
InputParameters
EpsKmodLubby2Templ<is_ad>::validParams()
{
  InputParameters params = RadialReturnStressUpdateTempl<is_ad>::validParams();
  params.addClassDescription("This class uses the stress update material in a radial return isotropic creep model"
    "This class computes the primary creep component of the modified Lubby2 creep.");

  // input parameters
  params.addRequiredParam<Real>("mvK", "Kelvin ViscoParameter");
  params.addRequiredParam<Real>("mk", "Kelvin Elastic Parameter");
  params.addRequiredParam<Real>("etaK0", "Initial Kelvin Viscosity");
  params.addRequiredParam<Real>("GK0", "Initial Kelvin Shear Modulus");
  return params;
}

template <bool is_ad>
EpsKmodLubby2Templ<is_ad>::EpsKmodLubby2Templ(
    const InputParameters & parameters)
  : RadialReturnStressUpdateTempl<is_ad>(parameters),
    _creep_strain(this->template declareGenericProperty<RankTwoTensor, is_ad>("creep_strain")),
    _creep_strain_old(this->template getMaterialPropertyOld<RankTwoTensor>("creep_strain")),
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
  _creep_strain[_qp].zero();
  _effective_creep_increment[_qp] = 0.0;
}

template <bool is_ad>
void
EpsKmodLubby2Templ<is_ad>::propagateQpStatefulProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
  _effective_creep_increment[_qp] = _effective_creep_increment_old[_qp];

  RadialReturnStressUpdateTempl<is_ad>::propagateQpStatefulPropertiesRadialReturn();
}

template <bool is_ad>
void
EpsKmodLubby2Templ<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & effective_trial_stress,
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
_effective_creep_increment[_qp] = _effective_creep_increment_old[_qp];
}

template <bool is_ad>
GenericReal<is_ad>
EpsKmodLubby2Templ<is_ad>::computeResidual(
    const GenericReal<is_ad> & effective_trial_stress, const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const GenericReal<is_ad> etaKMax = _etaK0 * std::exp(_mvK * stress_delta);
  const GenericReal<is_ad> GKMax = _GK0 * std::exp(_mk * stress_delta);
  _effective_creep_increment[_qp] = _effective_creep_increment_old[_qp] + scalar;
  const GenericReal<is_ad> creep_rate = stress_delta / (3.0 * etaKMax) - (GKMax*_effective_creep_increment[_qp])/(etaKMax);
  return creep_rate * _dt - scalar;
}

template <bool is_ad>
GenericReal<is_ad>
EpsKmodLubby2Templ<is_ad>::computeDerivative(
    const GenericReal<is_ad> & effective_trial_stress, const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  const GenericReal<is_ad> etaKMax = _etaK0 * std::exp(_mvK * stress_delta);
  const GenericReal<is_ad> GKMax = _GK0 * std::exp(_mk * stress_delta);
  _effective_creep_increment[_qp] = _effective_creep_increment_old[_qp] + scalar;
  const GenericReal<is_ad> creep_rate_derivative = ((_three_shear_modulus * GKMax *_effective_creep_increment[_qp]) / ( etaKMax)) * (_mvK - _mk);
  return creep_rate_derivative * _dt - 1.0;
}


template <bool is_ad>
void
EpsKmodLubby2Templ<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & creep_strain_increment)
{
  _creep_strain[_qp] += creep_strain_increment;
}


template class EpsKmodLubby2Templ<false>;
template class EpsKmodLubby2Templ<true>;
