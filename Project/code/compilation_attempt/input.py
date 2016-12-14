#####################################################################################################
# UNREFLECTED (BARE) SPHERE
#####################################################################################################
reflector=False
#####################################################################################################
# specify the radius of the core
core_radius = 10.0
# specify the number of cells
nCells = 20
# specify quadrature degree for mu and w
degree = 2
# total xs
sigma_tot = 2.5
# scattering xs
sigma_s = 0.5
# order of P_N expansion of scattering
L = 4
# incoming flux from boundary
hf=0.0
#---------------------------------------------------------------------------------------------------#
# external source: you can choose uniform or point source
sourcetype = 'uniform'
# sourcetype = 'point source'
# specify the value of that source
Q = 0.5
# if sourcetype == 'point source', you must specify the location of that point source.
# that location must be within the sphere
# pt_src_loc = 5
#####################################################################################################
# END OF THE INPUT
#####################################################################################################
