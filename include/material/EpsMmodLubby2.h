//* AJ

#pragma once

#include "RadialReturnStressUpdate.h"

/**
* This class uses the stress update material in a radial return isotropic creep
* model.  This class computes the primary creep component of the modified Lubby2 creep.
 */

template <bool is_ad>
class EpsMmodLubby2Templ : public RadialReturnStressUpdateTempl<is_ad>
{
public:
  static InputParameters validParams();

  EpsMmodLubby2Templ(const InputParameters & parameters);

  using Material::_qp;
  using Material::_dt;
  using RadialReturnStressUpdateTempl<is_ad>::_base_name;
  using RadialReturnStressUpdateTempl<is_ad>::_three_shear_modulus;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() override;

  virtual void
  computeStressInitialize(const GenericReal<is_ad> & effective_trial_stress,
                          const GenericRankFourTensor<is_ad> & elasticity_tensor) override;
  virtual GenericReal<is_ad> computeResidual(const GenericReal<is_ad> & effective_trial_stress,
                                             const GenericReal<is_ad> & scalar) override;
  virtual GenericReal<is_ad> computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                               const GenericReal<is_ad> & scalar) override;
  virtual void
  computeStressFinalize(const GenericRankTwoTensor<is_ad> & plastic_strain_increment) override;


  /// plastic strain in this model
  GenericMaterialProperty<RankTwoTensor, is_ad> & _plastic_strain;

  /// old value of plastic strain
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;

  /// Maxwell viscosity parameter
  const Real _mvM;

  /// Maxwell initial viscosity
  const Real _etaM0;

};

typedef EpsMmodLubby2Templ<false> EpsMmodLubby2;
typedef EpsMmodLubby2Templ<true> ADEpsMmodLubby2;
