cdef extern from "clf_mcmc.h":
	
	void getlf(double *param, double *phi_analysis)
	void getwp(double *param, double *w_analysis)

import numpy as np

def phi(para,analysisphi):

	if not para.flags['C_CONTIGUOUS']:
		para = np.ascontiguousarray(para)
	if not analysisphi.flags['C_CONTIGUOUS']:
		analysisphi = np.ascontiguousarray(analysisphi) 
   
	cdef double[::1] para_memview = para
	cdef double[::1] analysisphi_memview = analysisphi

	getlf(&para_memview[0], &analysisphi_memview[0])

	return analysisphi

def wp(para,analysiswp):

	if not para.flags['C_CONTIGUOUS']:
		para = np.ascontiguousarray(para)
	if not analysiswp.flags['C_CONTIGUOUS']:
		analysiswp = np.ascontiguousarray(analysiswp)
   
	cdef double[::1] para_memview = para
	cdef double[::1] analysiswp_memview = analysiswp

	getwp(&para_memview[0], &analysiswp_memview[0])

	return analysiswp

#from libc.stdio cimport printf