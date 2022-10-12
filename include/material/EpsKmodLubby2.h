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

/**
 * This class uses the stress update material in a radial return isotropic creep
 * model.  This class computes the secondary creep component of the modified Lubby2 creep.
 */


template <bool is_ad>
class EpsKmodLubby2Templ : public RadialReturnCreepStressUpdateBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  EpsKmodLubby2Templ(const InputParameters & parameters);

  virtual bool substeppingCapabilityEnabled() override;

  virtual void resetIncrementalMaterialProperties() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() override;

  virtual void
  computeStressInitialize(const GenericReal<is_ad> & effective_trial_stress,
                          const GenericRankFourTensor<is_ad> & elasticity_tensor) override{};  

virtual void
  computeStressFinalize(const GenericRankTwoTensor<is_ad> & plastic_strain_increment) override;


  virtual GenericReal<is_ad> computeResidual(const GenericReal<is_ad> & effective_trial_stress,
                                             const GenericReal<is_ad> & scalar) override
  {
    return computeResidualInternal<GenericReal<is_ad>>(effective_trial_stress, scalar);
  }


  virtual GenericReal<is_ad> computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                               const GenericReal<is_ad> & scalar) override;


  virtual GenericChainedReal<is_ad>
  computeResidualAndDerivative(const GenericReal<is_ad> & effective_trial_stress,
                               const GenericChainedReal<is_ad> & scalar) override
  {
    return computeResidualInternal<GenericChainedReal<is_ad>>(effective_trial_stress, scalar);
  }


    /// scalar effective strain rate in this model
    GenericMaterialProperty<Real, is_ad> & _effective_creep_increment;

    /// old value of scalar effective strain rate in this model
    const MaterialProperty<Real> & _effective_creep_increment_old;

    /// Kelvin ViscoParameter
    const Real _mvK;

    /// Kelvin Elastic Parameter
    const Real _mk;

    /// Initial Kelvin Viscosity
    const Real _etaK0;

    /// Initial Kelvin Shear Modulus
    const Real _GK0;

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

typedef EpsKmodLubby2Templ<false> EpsKmodLubby2;
typedef EpsKmodLubby2Templ<true> ADEpsKmodLubby2;
