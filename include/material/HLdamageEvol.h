#pragma once

#include "RadialReturnCreepStressUpdateBase.h"

/*
* This class takes advantage of the newton iteration in RadialReturnStressUpdate
to solve for the incremental value of a damage variable. Based on this incremental
value, the evolution of the damage variable is then solved iteratively.
*/


template <bool is_ad>
class HLdamageEvolTempl : public RadialReturnCreepStressUpdateBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  HLdamageEvolTempl(const InputParameters & parameters);

  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() override;


protected:
  virtual void
computeStressInitialize(const GenericReal<is_ad> & /*effective_trial_stress*/,
                          const GenericRankFourTensor<is_ad> & /*elasticity_tensor*/,
                          RankTwoTensor & stress_new,
                          const RankTwoTensor deviatoric_stress);

//  virtual void
//  computedamagestate(const GenericReal<is_ad> & scalar);

  //virtual void computeQpProperties() override;

  virtual void iterationFinalize(GenericReal<is_ad> scalar) override;

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

// model params
  Real _a4;
  Real _a5;
  Real _a6;
  Real _a7;
  Real _a8;
  Real _a9;
  Real _a10;
  Real _a15;
  Real _a17;

// 3rd deviatoric stress invariant and the min. principal stress
  Real J3_sig;
  Real smin;
  Real _scalar;

  // old and new damage property/variable being solved for:
  GenericMaterialProperty<Real, is_ad> & damage_variable;
  const MaterialProperty<Real> & damage_variable_old;

  using Material::computeQpProperties;
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

typedef HLdamageEvolTempl<false> HLdamageEvol;
typedef HLdamageEvolTempl<true> ADHLdamageEvol;
