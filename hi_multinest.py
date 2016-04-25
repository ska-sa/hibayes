#!/usr/bin/env python

"""
# hi_multinest.py

Bayesian Monte Carlo model fitting routine for global 21-cm cosmic
dawn data.

Jonathan Zwart
Danny Price
Gianni Bernardi
April 2016

Usage:

./hi_plot.py config_file.ini

"""


import os
import sys
import time
import numpy
from numpy import exp, log, sqrt, pi
import pymultinest
from scipy.special import erf
from hibayes.priors import Priors
from hibayes.profile_support import profile
from hibayes.constants import *
from hibayes.spectral_models import T_fg, T_HI, sigma
from hibayes.sky_model import generate_simulated_data, nu_1
from hibayes.parse_config import parse_config
from hibayes.utils import *

import pprint

# Turn off divide by zero errors
numpy.seterr(divide='ignore', invalid='ignore')

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print 'usage:'
        print
        print 'with MPI [NB mpiEXEC on marchmain]:'
        print '       mpiexec -n NPROC ./hi_multinest.py config_filename.ini'
        print
        sys.exit(0)

param_file = sys.argv[-1]
# load runtime parameters
rp = parse_config(param_file)
print "Runtime parameters"
pprint.pprint(rp)
time.sleep(2)


#-------------------------------------------------------------------------------

pri = Priors()

@profile
def myprior(cube, ndim, nparams):
    """
    """

    cube[0] = pri.UniformPrior(cube[0], -600.0, 0.0)  # A/mK
    cube[1] = pri.UniformPrior(cube[1], 40.0, 90.0)  # nu_HI/MHz
    cube[2] = pri.UniformPrior(cube[2], 0.0, 10.0)  # sigma_HI/MHz

    return

#-------------------------------------------------------------------------------

@profile
def myloglike(cube, ndim, nparams):
    """
    Form of lhood assumes uncertainties are gaussian and uncorrelated
    The latter can be extended later
    """

    A_HI     = cube[0]
    nu_HI    = cube[1]
    sigma_HI = cube[2]
    #c=[cube[ic+3] for ic in range(ndim)]

    loglike = 0.0
    for idatum in range(len(freqs)):
        #model = cube[0]*x[idatum]+cube[1]
        Tsky = T_HI(A_HI, nu_HI, sigma_HI, freqs[idatum], norm=True) + T_fg(nu_1, c, nc, freqs[idatum])
        sig = 1000 * sigma(Tsky, BW, tObs)
        #print idatum,freqs[idatum],Tsky,sig
        chisq = 0.5 * ((Tmeas[idatum] - Tsky) / sig) ** 2.0
        prefactor = 0.5 * log(2.0 * pi * sig ** 2.0)
        loglike -= prefactor + chisq
    #sys.exit(0)
    return loglike


#-------------------------------------------------------------------------------

#pri=Priors()
@profile
def mypriorj(cube, ndim, nparams, fg_only=False, bg_only=False, log_prior=False):
    """
    Joint prior - HI and bandpass
    """

    if log_prior:
        cube[0] = pri.GeneralPrior(cube[0], rp["A_HI_PRIOR"], rp["A_HI_MIN"], rp["A_HI_MAX"])  # comment for log priors
    else:
        cube[0] = pri.GeneralPrior(cube[0], 'U', rp["A_HI_MIN"], rp["A_HI_MAX"])  # A/mK  comment in for "normal runs"

    if not fg_only and not bg_only:
        cube[1] = pri.GeneralPrior(cube[1], 'U', rp["NU_MIN"], rp["NU_MAX"])
    elif not fg_only and bg_only:
        # GB: comment this out if you're fitting for the HI too
        cube[1]=pri.UniformPrior(cube[1], rp["FREQ_MIN"],rp["FREQ_MAX"]) # nu_HI/MHz
    elif fg_only and not bg_only:
        # GB: comment this out if you're fitting for the foregrounds only
        cube[1]=pri.GeneralPrior(cube[1],'U', rp["FREQ_MIN"], rp["FREQ_MAX"])

    cube[2] = pri.GeneralPrior(cube[2], 'U', rp["SIGMA_HI_MIN"], rp["SIGMA_HI_MAX"])  # sigma_HI/MHz
    #cube[2]=pri.UniformPrior(cube[2],0.01,0.01) # sigma_HI/MHz
    #cube[2]=pri.UniformPrior(cube[2],SIGMA_HI_TRUE,SIGMA_HI_TRUE) # sigma_HI/MHz
    #for i in range(ndim-3):
    #    #cube[i+3]=pri.UniformPrior(cube[i+3],-BP_PRIOR_RANGE,+BP_PRIOR_RANGE)
    #    cube[i+3]=pri.UniformPrior(cube[i+3],c[i],c[i])
    #cube[3]=pri.UniformPrior(cube[3],-BP_PRIOR_RANGE,+BP_PRIOR_RANGE)
    #cube[3]=pri.UniformPrior(cube[3],c[0],c[0])
    #cube[3]=pri.UniformPrior(cube[3],-BP_PRIOR_RANGE,+BP_PRIOR_RANGE)

    for ic in range(ndim - 3):
        cube[ic + 3] = pri.GeneralPrior(cube[ic + 3], 'U', -rp["BP_PRIOR_RANGE"], +rp["BP_PRIOR_RANGE"])
    
    #print "%2.2f | %2.2f | %2.2f" % (cube[0], cube[1], cube[2])
    return


#-------------------------------------------------------------------------------

@profile
def myloglikej(cube, ndim, nparams):
    """
    Form of lhood assumes uncertainties are gaussian and uncorrelated
    The latter can be extended later
    
    """
    #numpy.savetxt('Tmeas.txt',Tmeas)	# GB: just a test to see that it reads the correct spectrum when run on real data
    #numpy.savetxt('freqs.txt',freqs)	# GB: just a test to see that it reads the correct frequency when run on real data

    A_HI = cube[0]
    nu_HI = cube[1]
    sigma_HI = cube[2]
    #c=[cube[ic+3] for ic in range(ndim)]
    #print A_HI

    c_fit = [cube[i + 3] for i in range(ndim - 3)]
    #print c

    loglike = 0.0
    for idatum in range(len(freqs)):
        Tsky = 1.0e-3 * T_HI(A_HI, nu_HI, sigma_HI, freqs[idatum]) + T_fg(nu_1, c_fit, rp["nc_fit"], freqs[idatum])
        #print idatum,freqs[idatum],Tsky,Tp
        #print idatum,Tmeas[idatum],Tf(nu_1,c,nc,freqs[idatum]),1.0e-3*T_HI(A_HI,nu_HI,sigma_HI,freqs[idatum]),Tsky,sig,c
        #sig=sigma(Tf(nu_1,c,nc,freqs[idatum]),BW,tObs)
        if rp["spectrum_errors"] is not None:
            sig = errors[idatum]
        else:
            sig = sigma(Tsky, rp["BW"], rp["tObs"])
            #print idatum,Tpure[idatum],Tsky
            #print freqs[idatum],Tpure[idatum],Tmeas[idatum],Tsky,Tpure[idatum+1]-Tsky,sig,Thi[idatum]#,c
        #        print Thi[idatum],1.0e-3*T_HI(A_HI,nu_HI,sigma_HI,freqs[idatum])
        chisq = 0.5 * ((Tmeas[idatum] - Tsky) / sig) ** 2.0
        prefactor = 0.5 * log(2.0 * pi * sig ** 2.0)
        loglike -= prefactor + chisq

    #sys.exit(0)

    return loglike


#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------

@profile
def main():
    """
    """

    #try:
    #    execfile(param_file)
    #except IOError:
    #    from settings import *

    global freqs, Tmeas, FREQ_MIN, FREQ_MAX, n_params

    if not os.path.exists(rp["outdir"]):
        try:
            os.mkdir(rp["outdir"])
        except:
            pass


    if rp["ledaSpec"] is None:
        Tmeas, freqs = generate_simulated_data(rp)
    else:
        freqs = numpy.genfromtxt(rp["ledaFreqs"])
        Tmeas = numpy.genfromtxt(rp["ledaSpec"])
        FREQ_MIN = freqs[0];
        FREQ_MAX = freqs[-1]

    n_params = rp["nc_fit"] + 3


    progress = pymultinest.ProgressPlotter(n_params=n_params,  interval_ms=10000,
                                           outputfiles_basename=rp["outputfiles_basename"])
    #progress.start()
    
    pymultinest.run(myloglikej, mypriorj, n_params, resume=False, verbose=True,
                    multimodal=rp["multimodal"], max_modes=rp["max_modes"], write_output=True,
                    n_live_points=rp["n_live_points"],
                    evidence_tolerance=rp["evidence_tolerance"],
                    mode_tolerance=rp["mode_tolerance"], seed=rp["seed"], max_iter=rp["max_iter"],
                    importance_nested_sampling=rp["do_ins"],
                    outputfiles_basename=rp["outputfiles_basename"], init_MPI=False)

    #progress.stop()

    #contour_plot.contourTri(pylab.loadtxt('chains_test/1-post_equal_weights.dat'),
    # line=True,outfile='chains_test/test.png',col=('red','blue'),labels=parameters,
    # binsize=50,truth=plotTruth,ranges=plotRanges,autoscale=True)

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    # Run Multinest
    ret = main()
    sys.exit(ret)






