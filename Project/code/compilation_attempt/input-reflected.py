#####################################################################################################
# REFLECTED (BARE) SPHERE
#####################################################################################################
# specify the radius of the core
core_radius = 5.0
# reflector thickness
refl_thickness = 10.0
# specify the number of cells
nCells_core = 100
nCells_refl = 100
# Angular variables - quadrature degree
degree = 2
# Neutronic variables - specify the XS in each region
# Core
sigma_tot_core = 2.0
sigma_s_core = 1.0
# Reflector
sigma_tot_refl = 5.0
sigma_s_refl = 2.0
# external source
Q = 0.8
# Flux incoming from outside the sphere
hf=0
# order of P_N expansion of scattering
L = 4
