#####################################################################################################
import input
import functions
from functions import initialize_inputs
from functions import iterate_flux
import os
import re
#####################################################################################################
outfile = open('output.txt', 'w')
outfile.write('Solution to Transport Equation in 1-D Spherical Coordinates \n')
outfile.write('By Milos Atz \n')
outfile.write('NE255 Final Project Fall 2016\n')
outfile.write('#####################################################################################################\n')
outfile.write('Inputs: \n')
inputs = open('input.py', 'r').read()
outfile.write(inputs)
outfile.write('\n')
outfile.write('#####################################################################################################\n')
#####################################################################################################
cc, cs, mu, w, hf, extsource, L, xs_tot, xs_s = initialize_inputs()
#try:
sf, af = iterate_flux(cc, cs, mu, w, hf, extsource, L, xs_tot, xs_s)
#print(functions.init_wdd(cc, cs, mu[0], extsource, xs_tot, hf))

#except:
#	print('not converged')
#outfile.write('Results : \n')
#outfile.write('cell centers (positions) \n')
#outfile.write(str(cc)+'\n')
#outfile.write('scalarflux \n')
#outfile.write(str(sf)+'\n')
#outfile.write('angularflux \n')
#for m in af:
#	outfile.write('mu = '+str(m)+'\n')
#	outfile.write('angular flux (mu) = \n')
#	outfile.write(af[m])