//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "RadialReturnCreepStressUpdateBase.h"
#include<algorithm>

/**
 * This class uses the stress update material in a radial return isotropic creep
 * model.
 */

template <bool is_ad>
class HouLuxTempl : public RadialReturnCreepStressUpdateBaseTempl<is_ad>
{
public:
  static InputParameters validParams();
  HouLuxTempl(const InputParameters & parameters);
  virtual bool substeppingCapabilityEnabled() override;
  virtual void resetIncrementalMaterialProperties() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() override;


  virtual void
  computeStressInitialize(const GenericReal<is_ad> & /*effective_trial_stress*/,
                          const GenericRankFourTensor<is_ad> & /*elasticity_tensor*/,
                          RankTwoTensor & stress_new,
                          const RankTwoTensor deviatoric_stress);
  virtual void
  computeStressFinalize(const GenericRankTwoTensor<is_ad> & plastic_strain_increment) override;

  //Declare residuals when automatic_differentiation = false

  virtual GenericReal<is_ad>
  computeResidual(const GenericReal<is_ad> & effective_trial_stress,
                                             const GenericReal<is_ad> & scalar) override
  {
    return computeResidualInternal<GenericReal<is_ad>>(effective_trial_stress, scalar);
  }

  /// Declare Derivatives when automatic_differentiation = false
  virtual GenericReal<is_ad>
  computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                               const GenericReal<is_ad> & scalar) override;

  /// Declare computeResidual And Derivatives for automatic_differentiation:
  virtual GenericChainedReal<is_ad>
  computeResidualAndDerivative(const GenericReal<is_ad> & effective_trial_stress,
                              const GenericChainedReal<is_ad> & scalar) override
  {
  return computeResidualInternal<GenericChainedReal<is_ad>>(effective_trial_stress, scalar);
  }



   /**
   * Damage parameter updated according to equation (20) in:
   * Mechanical and hydraulic behavior of rock salt in the excavation
   * disturbed zone around underground facilities.
   * Hou (2003)
   */

   //virtual GenericReal<is_ad>
   //updateDamageParam (const GenericReal<is_ad> & effective_trial_stress,
    //                                         const GenericReal<is_ad> & scalar);

  /// Declare damage parameters: new and old
    GenericMaterialProperty<Real, is_ad> & _damage_param;
    const MaterialProperty<Real> & _damage_param_old;

  /// Maxwell initial viscosity
  const Real _etaM0;
  /// Maxwell viscosity parameter
  const Real _mvM;
  /// Kelvin ViscoParameter
  const Real _mvK;
  /// Kelvin Elastic Parameter
  const Real _mk;
  /// Initial Kelvin Viscosity
  const Real _etaK0;
  /// Initial Kelvin Shear Modulus
  const Real _GK0;

  /// Declare scalar effective kelvin strain rate
  GenericMaterialProperty<Real, is_ad> & _kelvin_creep_rate;
  const MaterialProperty<Real> & _kelvin_creep_rate_old;


  // model params
    Real _a1;
    Real _a2;
    Real _a3;
    Real _a4;
    Real _a5;
    Real _a6;
    Real _a7;
    Real _a8;
    Real _a9;
    Real _a10;
    Real _a11;
    Real _a12;
    Real _a13;
    Real _a15;
    Real _a16;
    Real _a17;
    Real _F_star;
    Real _evol;

  // 3rd deviatoric stress invariant and the min. principal stress
    Real J3_sig;
    Real smin;


  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_qp;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_dt;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_t;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_three_shear_modulus;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_creep_strain;
  using RadialReturnCreepStressUpdateBaseTempl<is_ad>::_creep_strain_old;

private:

template <typename ScalarType>
  ScalarType computeResidualInternal(const GenericReal<is_ad> & effective_trial_stress,
                                     const ScalarType & scalar);

};

typedef HouLuxTempl<false> HouLux;
typedef HouLuxTempl<true> ADHouLux;
