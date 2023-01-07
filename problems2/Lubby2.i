# Heusserman et. al.

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
    type = ADMaterialRealAux
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
    use_automatic_differentiation = true
  []
[]

[Functions]
  [top_pull]
   type = PiecewiseLinear
   x = ' 0  35 ' #  should be x= '0 1 15' where loading=0 to 1 and creep-phase=1 to 15.   
   y = ' 0 -12e6' #  should be y= '0 -2.5e5 -2.5e5', where 0 to -2.5e5 indicates loading phase, -2.5e5 to -2.5e5 indicates creep phase where load is kept constant. 
  []
[]

[BCs]             # if DirichletBC given a value equal 0, it means fixed (pin). If no value, it means roller. ideally, for uniaxial compressive, all sides should be roller,except fixed bottom and top pull.
  [u_top_pull]
    type = Pressure
    variable = disp_y
    boundary = top
    factor = 1                     #-2.5e5
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
    type = ADComputeIsotropicElasticityTensor
    block = 0
    youngs_modulus = 2e6 # 1.58e9 # 2e6 #   use e4
    poissons_ratio = 0.3 #0.335  # 0.3 #
  []
  [radial_return_stress]
    type = ADComputeMultipleInelasticStress
     block = 0
    inelastic_models = "_L2creep"
    combined_inelastic_strain_weights = '1.0'
  []
  [_L2creep]
    type = ADLubby2
    block = 0
    etaM0 = 1.21e8          #4.03e7     #
    mvM   = -0.327          # -0.327E-6 #
    GK0   = 1.88e5          #6.27e4     #
    mk    = -0.254          # -0.254e-6 #
    etaK0 = 4.98e5          #1.66e5     #
    mvK   = -0.267          #-0.267e-6  #
    EQstress = "12 0 0  0 0  0  0 0 0"
    internal_solve_full_iteration_history = true
#    use_substep = false
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK' 

  petsc_options = '-snes_ksp'
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = '101'
  line_search = 'none'
#  automatic_scaling = false
# computing_scaling_once = false
  l_max_its = 20
  nl_max_its = 20
  nl_rel_tol = 1e-8  #1e-6
  nl_abs_tol = 1e-11 #1e-6
  l_tol = 1e-5       #1e-5
  start_time = 0.0
  end_time = 35.0

#  [TimeIntegrator]
#  type = CrankNicolson
#  []
 # num_steps = 35
 # dt = 1

 # [./TimeStepper]
 #   type = FunctionDT
 #   function = top_pull
 # [../]
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
 # exodus = true
  [csv]
  type = CSV
  execute_on = 'initial timestep_end'
  []
[]

#[Debug]
#  show_material_props = true
#  internal_solve_full_iteration_history = true
#[]