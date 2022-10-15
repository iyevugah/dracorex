# 1x1x1 unit cube with uniform pressure on top face

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
[]

[AuxVariables]
  [damage_param]
    order = CONSTANT
    family = MONOMIAL
  []
   [kelvin_creep_rate]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
   [damage_evol]
    type = MaterialRealAux
    variable = damage_param
    property = damage_param
    execute_on = 'initial timestep_end'
  []
    [kelvin_creep_rate_evol]
    type = MaterialRealAux
    variable = kelvin_creep_rate
    property = kelvin_creep_rate
    execute_on = 'initial timestep_end'
  []
[]


[Modules/TensorMechanics/Master]
  [all]
    strain = SMALL
    incremental = true
    add_variables = true 
    generate_output = 'stress_yy creep_strain_xx creep_strain_yy creep_strain_zz elastic_strain_yy' 
  []
[]

[Functions]
  [top_pull]
    type = PiecewiseLinear
    x = '0 1'
    y = '1 1'
  []
[]

[BCs]
  [u_top_pull]
    type = Pressure
    variable = disp_y
    boundary = top
    factor = '-10.0e6'
    function = top_pull
  []
  [u_bottom_fix]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  []
  [u_yz_fix]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  []
  [u_xy_fix]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  []
[]


[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2e10     # 2e11
    poissons_ratio = 0.3      # 0.3
  []
  [radial_return_stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'HouLux_Eps'
    tangent_operator = elastic
    combined_inelastic_strain_weights = '1.0'
  []
  [HouLux_Eps]
    type = HouLuxEps
    mvM =  -2.47e-6             #-2.67e-8  1.9e-6     I scaled the model parameter by e-7 
    etaM0 = 2.03e7           # 4e7      2.03e7  
    mvK =   -1.68e-6            # -3.27e-8
    mk =    -1.91e-6           # -2.54e-8
    etaK0 =  8.94e3            # 1.66e5
    GK0 =    5.08e5          # 6.27e4
      a4 = 0.3
      a5 = 0.05
      a6 = 67.0
      a7 = 41.0
      a8 = 0.25
      a9 = 1.0
     a10 = 0.25
     a15 = 1.67e-1
     a16 = 1.0e-8
     a17 = 3.5e-8
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK' 

  petsc_options = '-snes_ksp'
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = '101'
  line_search = 'none'
  automatic_scaling = true
  computing_scaling_once = false

  l_max_its = 20
  nl_max_its = 20
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_tol = 1e-5
  start_time = 0.0
  end_time = 100.0
  num_steps = 100
  dt = 0.001
[]

[Postprocessors]
    [elastic]
    type = PointValue
    point = '1 1 1'
    variable = elastic_strain_yy
  []
    [creep_strainXX]
    type = PointValue
    point = '1 1 1'
    variable = creep_strain_xx
  []
    [creep_strainYY]
    type = PointValue
    point = '1 1 1'
    variable = creep_strain_yy
  []
    [creep_strainZZ]
    type = PointValue
    point = '1 1 1'
    variable = creep_strain_zz
  []
    [stressYY]
    type = PointValue
    point = '1 1 1'
    variable = stress_yy
  []
  [damage_evol]
   type = PointValue
   point = '1 1 1'
   variable = damage_param
  []
  [scalar_Kelvin_strain_rate]
   type = PointValue
   point = '1 1 1'
   variable = kelvin_creep_rate
  []
[]

[Outputs]
  exodus = true
  [csv]
  type = CSV
  execute_on = 'initial timestep_end'
  []
[]