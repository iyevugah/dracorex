//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HouLux.h"

registerMooseObject("dracorexApp", HouLux);
registerMooseObject("dracorexApp", ADHouLux);

template <bool is_ad>
InputParameters
HouLuxTempl<is_ad>::validParams()
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
  params.addRequiredParam<Real>("a1", "model parameter 1");
  params.addRequiredParam<Real>("a2", "model parameter 2");
  params.addRequiredParam<Real>("a3", "model parameter 3");
  params.addRequiredParam<Real>("a4", "model parameter 4");
  params.addRequiredParam<Real>("a5", "model parameter 5");
  params.addRequiredParam<Real>("a6", "model parameter 6");
  params.addRequiredParam<Real>("a7", "model parameter 7");
  params.addRequiredParam<Real>("a8", "model parameter 8");
  params.addRequiredParam<Real>("a9", "model parameter 9");
  params.addRequiredParam<Real>("a10", "model parameter 10");
  params.addRequiredParam<Real>("a11", "model parameter 11");
  params.addRequiredParam<Real>("a12", "model parameter 12");
  params.addRequiredParam<Real>("a13", "model parameter 13");
  params.addRequiredParam<Real>("a15", "model parameter 15");
  params.addRequiredParam<Real>("a16", "model parameter 16");
  params.addRequiredParam<Real>("a17", "model parameter 17");
  params.addRequiredParam<Real>("F_star", "model parameter F_star");
  params.addRequiredParam<Real>("evol", "model parameter evolution");
  return params;
}

template <bool is_ad>
HouLuxTempl<is_ad>::HouLuxTempl(const InputParameters & parameters)
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
    _a1 (this->template getParam<Real>("a1")),
    _a2 (this->template getParam<Real>("a2")),
    _a3 (this->template getParam<Real>("a3")),
    _a4 (this->template getParam<Real>("a4")),
    _a5 (this->template getParam<Real>("a5")),
    _a6 (this->template getParam<Real>("a6")),
    _a7 (this->template getParam<Real>("a7")),
    _a8 (this->template getParam<Real>("a8")),
    _a9 (this->template getParam<Real>("a9")),
    _a10 (this->template getParam<Real>("a10")),
    _a11 (this->template getParam<Real>("a11")),
    _a12 (this->template getParam<Real>("a12")),
    _a13 (this->template getParam<Real>("a13")),
    _a15 (this->template getParam<Real>("a15")),
    _a16 (this->template getParam<Real>("a16")),
    _a17 (this->template getParam<Real>("a17")),
    _F_star (this->template getParam<Real>("F_star")),
    _evol (this->template getParam<Real>("evol"))

{
  if (_etaM0 == 0.0 && _etaK0 == 0.0 &&  _a3 == 0 &&  _evol == 0)
    mooseError("HouLux: at least one of the creep should be active.");
}


template <bool is_ad>
void
HouLuxTempl<is_ad>::initQpStatefulProperties()
{
  _kelvin_creep_rate[_qp] = 0.0;
  _damage_param[_qp] = 0.0; // Initialize damage parameter
}

template <bool is_ad>
void
HouLuxTempl<is_ad>::propagateQpStatefulProperties()
{
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp];
  _damage_param[_qp] = _damage_param_old[_qp];
  RadialReturnStressUpdateTempl<is_ad>::propagateQpStatefulPropertiesRadialReturn();
}


template <bool is_ad>
void
HouLuxTempl<is_ad>::computeStressInitialize(const GenericReal<is_ad> & /*effective_trial_stress*/,
                                               const GenericRankFourTensor<is_ad> & /*elasticity_tensor*/,
                                          RankTwoTensor & stress_new, const RankTwoTensor deviatoric_stress)
{
  _kelvin_creep_rate[_qp] = _kelvin_creep_rate_old[_qp];
  _damage_param[_qp] = _damage_param_old[_qp];
  Real J3_sig = deviatoric_stress.det(); // 3rd invariant of deviatoric stress
  Real smin = ((stress_new) - deviatoric_stress).trace(); // minimum principal stress
}


/// Compute the Residuals whether automatic_differentiation = false or true
/// in the SingleVariableReturnMappingSolution.C file in MOOSE
template <bool is_ad>
template <typename ScalarType>
ScalarType
HouLuxTempl<is_ad>::computeResidualInternal(const GenericReal<is_ad> & effective_trial_stress,
                                                 const ScalarType & scalar)
{
  const ScalarType stress_delta = effective_trial_stress - (_three_shear_modulus * scalar);
  //computing the Damage Parameter
  ScalarType lode_param = (27.0/2.0) * (J3_sig) / pow(stress_delta,3.0);  // Lode Parameter
  ScalarType lode_ang = 1.0 - ((2.0/pi) * acos(lode_param));       // Lode Angle
  ScalarType eta_D = (1.0) - (_a4 * exp(-_a5*smin));
  ScalarType beta_TC = (_a6) - (_a7 * exp(-_a8*smin));
  ScalarType kappa_beta = pow(1.0 / (cos(lode_ang + (pi / 6.0)) + _a9 * sin(lode_ang + (pi / 6.0))),
                                exp(-(_a10) * smin));
  ScalarType F_ds = stress_delta - ((eta_D)*(beta_TC)*(kappa_beta));
  ScalarType mac_Fds = std::max((F_ds / _F_star), 1.0);
  ScalarType F_dz = 6.0 * (-smin);
  ScalarType mac_Fdz = std::max((F_dz / _F_star), 1.0);
  ScalarType F_h = ((2.0*smin)/3.0) + (2.0/(3.0*_a5)*log((_a6-stress_delta)/(_a6)));
  ScalarType mac_Fh = std::max((F_h / _F_star), 1.0);
  _damage_param[_qp] = _damage_param_old[_qp];
  //const ScalarType  damage_rate = (_a15/pow((1.0 - _damage_param[_qp]),_a17)) * pow((mac_Fds + mac_Fdz),_a16);
  const ScalarType  damage_rate = (_a15/(pow(1.0 - _damage_param[_qp],_a17 ))) * pow((F_ds + F_dz),_a16);
  _damage_param[_qp] = _damage_param_old[_qp] + MetaPhysicL::raw_value(damage_rate) * _dt; // update damage param

  const ScalarType etaM = _etaM0 * std::exp(_mvM * stress_delta);
  const ScalarType etaK = _etaK0 * std::exp(_mvK * stress_delta);
  const ScalarType GK = _GK0 * std::exp(_mk * stress_delta);
  _kelvin_creep_rate[_qp] =  _kelvin_creep_rate_old[_qp] + MetaPhysicL::raw_value(scalar);

  if (_etaM0 != 0.0 && _etaK0 != 0.0 && _a3 != 0 &&  _evol != 0)
  {
    // Maxwell, Kelvin, damage and healing residuals
  const ScalarType M_creep_rate = (stress_delta / ((1.0 -_damage_param[_qp]) * etaM));
  const ScalarType K_creep_rate = (stress_delta/((1.0 -_damage_param[_qp])* etaK)) -
                                   ((GK*(_kelvin_creep_rate[_qp]))/ ((std::sqrt(2.0/3.0))* etaK));
  const ScalarType D_creep_rate = (_a3/pow((1.0 -_damage_param[_qp]),_a2)) * (pow(mac_Fds,_a1) + pow(mac_Fdz,_a1));
  const ScalarType H_creep_rate = _evol * ((1/ (_a11 + _a12 * exp(_a13 *_evol))) * mac_Fh);
return (M_creep_rate + K_creep_rate + D_creep_rate + H_creep_rate) * _dt - scalar;

   }
       else if (_etaM0 != 0.0 &&  _etaK0 == 0.0 && _a3 == 0 &&  _evol == 0)
      // Maxwell
        {
          const ScalarType creep_rate = (stress_delta / ((1.0 -_damage_param[_qp]) * etaM));
         return creep_rate * _dt - scalar;
        }
        else if (_etaM0 == 0.0 &&  _etaK0 != 0.0 && _a3 == 0 &&  _evol == 0)
        // Kelvin
         {
          const ScalarType creep_rate = (stress_delta/((1.0-_damage_param[_qp])* etaK)) -
                                       ((GK*(_kelvin_creep_rate[_qp]))/ ((std::sqrt(2.0/3.0))* etaK));
         return creep_rate * _dt - scalar;
        }
       else if (_etaM0 == 0.0 &&  _etaK0 == 0.0 && _a3 != 0 &&  _evol == 0)
       // Damage
        {
        const ScalarType creep_rate  = (_a3/pow((1.0 -_damage_param[_qp]),_a2)) * (pow(mac_Fds,_a1) + pow(mac_Fdz,_a1));
        return creep_rate * _dt - scalar;
        }
        else if (_etaM0 == 0.0 &&  _etaK0 == 0.0 && _a3 != 0 &&  _evol == 0)
         // Maxwell & Kelvin
         {
           const ScalarType M_creep_rate = (stress_delta / ((1.0 -_damage_param[_qp]) * etaM));
           const ScalarType K_creep_rate = (stress_delta/((1.0-_damage_param[_qp])* etaK)) -
                                            ((GK*(_kelvin_creep_rate[_qp]))/ ((std::sqrt(2.0/3.0))* etaK));
          return (M_creep_rate + K_creep_rate) * _dt - scalar;
         }
         else if (_etaM0 == 0.0 &&  _etaK0 == 0.0 && _a3 != 0 &&  _evol == 0)
         // Maxwell, Kelvin & Damage
         {
           const ScalarType M_creep_rate = (stress_delta / ((1.0 -_damage_param[_qp]) * etaM));
           const ScalarType K_creep_rate = (stress_delta/((1.0-_damage_param[_qp])* etaK)) -
                                            ((GK*(_kelvin_creep_rate[_qp]))/ ((std::sqrt(2.0/3.0))* etaK));
           const ScalarType D_creep_rate = (_a3/pow((1.0 -_damage_param[_qp]),_a2)) * (pow(mac_Fds,_a1) + pow(mac_Fdz,_a1));
        return (M_creep_rate + K_creep_rate + D_creep_rate) * _dt - scalar;
      }
        else if (_etaM0 == 0.0 &&  _etaK0 == 0.0 && _a3 != 0 &&  _evol == 0)
        //  Damage & Heal
        {
          const ScalarType D_creep_rate = (_a3/pow((1.0 -_damage_param[_qp]),_a2)) * (pow(mac_Fds,_a1) + pow(mac_Fdz,_a1));
            const ScalarType H_creep_rate = _evol * ((1/ (_a11 + _a12 * exp(_a13 *_evol))) * (F_h/_F_star));
        return ( D_creep_rate + H_creep_rate) * _dt - scalar;
      }
  mooseError("HouLux: The total creep rate combination is unacceptable.");
 }


 // No need to compute derivatives. Set automatic_differentiation = true
 // in the SingleVariableReturnMappingSolution.C file in MOOSE
 template <bool is_ad>
 GenericReal<is_ad>
 HouLuxTempl<is_ad>::computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                             const GenericReal<is_ad> & scalar)
    {
         return 0.0;
    }


template <bool is_ad>
void
HouLuxTempl<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _creep_strain[_qp] += plastic_strain_increment;
}

template <bool is_ad>
void
HouLuxTempl<is_ad>::resetIncrementalMaterialProperties()
{
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <bool is_ad>
bool
HouLuxTempl<is_ad>::substeppingCapabilityEnabled()
{
  return this->template getParam<bool>("use_substep");
}

template class HouLuxTempl<false>;
template class HouLuxTempl<true>;

template Real HouLuxTempl<false>::computeResidualInternal<Real>(const Real &, const Real &);
template ADReal HouLuxTempl<true>::computeResidualInternal<ADReal>(const ADReal &,
                                                                     const ADReal &);
template ChainedReal
HouLuxTempl<false>::computeResidualInternal<ChainedReal>(const Real &, const ChainedReal &);

template ChainedADReal
HouLuxTempl<true>::computeResidualInternal<ChainedADReal>(const ADReal &, const ChainedADReal &);
