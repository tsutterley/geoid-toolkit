#!/usr/bin/env python
u"""
gauss_weights.py
Adapted by Tyler Sutterley (05/2013)

Computes the Gaussian weights as a function of degree
A normalized version of Jekeli's Gaussian averaging function

CALLING SEQUENCE:
	wl = gauss_weights(hw, lmax)

INPUTS:
	hw: Gaussian smoothing radius in kilometers
		Radius r corresponds to the distance at which the weight
		drops to half its peak value at the shortest wavelength
	lmax: Maximum degree of Stokes coefficients

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python (http://www.numpy.org)

REFERENCE:
	Jekeli, "Alternative Methods to Smooth the Earth's Gravity Field", 1981
		http://www.geology.osu.edu/~jekeli.1/OSUReports/reports/report_327.pdf

	Wahr, Molenaar and Frank, "Time variability of the Earth's gravity field:
		Hydrological and oceanic effects and their possible detection using
		GRACE", Journal of Geophysical Research: Solid Earth, 103(B12),
		30205-30229, http://dx.doi.org/10.1029/98JB02844

UPDATE HISTORY:
	Written 05/2013
"""
import numpy as np

def gauss_weights(hw, lmax):
	#-- allocate for output weights
	wl=np.zeros((lmax+1))
	#-- radius of the Earth in km
	rad_e=6371.0
	if (hw < 1.0e-10):
		#-- distance is smaller than cutoff
		wl[:]=1.0/(2.0*np.pi)
	else:
		#-- calculate gaussian weights using recursion
		b = np.log(2.0)/(1.0 - np.cos(hw/rad_e))
		#-- weight for degree 0
		wl[0] = 1.0/(2.0*np.pi)
		#-- weight for degree 1
		wl[1] = wl[0]*((1.0 +np.exp(-2.0*b))/(1. -np.exp(-2.0*b))-1.0/b)
		#-- valid flag
		valid = True
		#-- spherical harmonic degree
		l = 2
		#-- while valid (within cutoff)
		#-- and spherical harmonic degree is less than lmax
		while (valid and (l <= lmax)):
			#-- calculate weight with recursion
			wl[l] = (1.0-2.0*l)/b*wl[l-1]+wl[l-2]
			#-- weight is less than cutoff
			if (np.abs(wl[l]) < 1.0e-10):
				#-- set all weights to cutoff
				wl[l:lmax+1] = 1.0e-10
				#-- set valid flag
				valid = False
			#-- add 1 to l
			l += 1
	#-- return the gaussian weights
	return wl
