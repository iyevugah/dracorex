//* AJ

#pragma once

#include "RadialReturnStressUpdate.h"

/**
* This class uses the stress update material in a radial return isotropic creep
* model.  This class computes the primary creep component of the modified Lubby2 creep.
 */

template <bool is_ad>
class EpsKmodLubby2Templ : public RadialReturnStressUpdateTempl<is_ad>
{
public:
  static InputParameters validParams();

  EpsKmodLubby2Templ(const InputParameters & parameters);

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
  computeStressFinalize(const GenericRankTwoTensor<is_ad> & creep_strain_increment) override;


  /// plastic strain in this model
  GenericMaterialProperty<RankTwoTensor, is_ad> & _creep_strain;

  /// old value of plastic strain
  const MaterialProperty<RankTwoTensor> & _creep_strain_old;

  /// plastic strain in this model
  GenericMaterialProperty<Real, is_ad> & _effective_creep_increment;

  /// old value of plastic strain
  const MaterialProperty<Real> & _effective_creep_increment_old;

  /// Kelvin ViscoParameter
  const Real _mvK;

  /// Kelvin Elastic Parameter
  const Real _mk;

  /// Initial Kelvin Viscosity
  const Real _etaK0;

  /// Initial Kelvin Shear Modulus
  const Real _GK0;

};

typedef EpsKmodLubby2Templ<false> EpsKmodLubby2;
typedef EpsKmodLubby2Templ<true> ADEpsKmodLubby2;
