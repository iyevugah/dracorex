# 3x3x3 unit cube with uniform pressure on top face

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 3
  ny = 3
  nz = 3
[]

[AuxVariables]
   [kelvin_creep_rate]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
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
    use_automatic_differentiation = false
  []
[]

[Functions]
  [top_pull]
    type = PiecewiseLinear
    x = '0 15'
    y = '1 1'
  []
[]

[BCs]
  [u_top_pull]
    type = Pressure
    variable = disp_y
    boundary = top
    factor = -2.0e5
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
    block = 0
    youngs_modulus = 2e8
    poissons_ratio = 0.3
  []
  [radial_return_stress]
    type = ComputeMultipleInelasticStress
     block = 0
    inelastic_models = "Esp_modLubby2"
    tangent_operator = elastic
    combined_inelastic_strain_weights = '1.0'
  []
  [Esp_modLubby2]
# The following material params were taken from Zhang. et al. (2003) 
#Error-controlled implicit time integration of
#elasto-visco-plastic constitutive models for rock salt
    type = modLubby2
    block = 0
    etaM0 =  4.03e7           
    mvM   = -0.327e-6    #scaled by e-7
    GK0   = 6.27e4    
    mk    = -0.254e-6   #scaled by e-7
    etaK0 = 1.66e5    
    mvK   =  -0.267e-6   #scaled by e-7
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
  end_time = 15.0
 # num_steps = 10
 # dt = 0.1
[]

[Postprocessors]
    [elastic]
    type = PointValue
    point = '1 1 1'
    variable = elastic_strain_yy
  []
    [stressYY]
    type = PointValue
    point = '1 1 1'
    variable = stress_yy
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