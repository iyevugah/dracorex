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
    x = '0 100'
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
    youngs_modulus = 2e11     # 2e11
    poissons_ratio = 0.3      # 0.3
  []
  [radial_return_stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'LuxWolters_Eps'
    tangent_operator = elastic
    combined_inelastic_strain_weights = '1.0'
  []
  [LuxWolters_Eps]
    type = LuxWoltersEps    
    etaM0 =  3.0e7           
    mvM   = -3.621e-8    #scaled by e-7
    GK0   = 3.6342e4    
    mk    = -0.1831e-8   #scaled by e-7
    etaK0 = 6.64768e5    
    mvK   =  -2.898e-8   #scaled by e-7       
       L  = 0.055
      L1  = 0.0
       T  = 200
       a  = 0
       b  = 0 
      a4  = 0.82
      a5  = 0.1
      a6  = 63.0
      a7  = 32.0
      a8  = 0.22
      a9 = 0.17 
     a10 = 1e-3 
     a15 = 34e-3 
     a16 = 3.1e-3        #scaled by e-3
     a17 = 1.0 
   sigma0 = 1.0
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
  start_time = '0.0'
  end_time = '100'
  num_steps = 100
  dt = 0.01
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