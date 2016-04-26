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

./hi_multinest.py config_file.ini

"""

import os
import sys
import time
import numpy
import pymultinest
from mpi4py import MPI

from hibayes.profile_support import profile
from hibayes.parse_config import parse_config
from hibayes.likelihoods import logprior,loglike

import pprint

# Turn off divide by zero errors
numpy.seterr(divide='ignore', invalid='ignore')

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print 'usage:'
        print
        print 'with MPI:'
        print '       mpiexec -n NPROCS ./hi_multinest.py config_filename.ini'
        print
        sys.exit(0)

param_file = sys.argv[-1]
# load runtime parameters
rp = parse_config(param_file)
print "Runtime parameters"
pprint.pprint(rp)
time.sleep(2)

#-------------------------------------------------------------------------------

@profile
def main():
    """
    """

    if not os.path.exists(rp["outdir"]):
        try:
            os.mkdir(rp["outdir"])
        except:
            pass

    n_params = rp["nc_fit"] + 3

    #progress = pymultinest.ProgressPlotter(n_params=n_params,  interval_ms=10000,
    #                                       outputfiles_basename=rp["outputfiles_basename"])
    #progress.start()
    
    pymultinest.run(loglike, logprior, n_params, resume=False, verbose=True,
                    multimodal=rp["multimodal"], max_modes=rp["max_modes"], write_output=True,
                    n_live_points=rp["n_live_points"],
                    evidence_tolerance=rp["evidence_tolerance"],
                    mode_tolerance=rp["mode_tolerance"],
                    seed=rp["seed"],
                    max_iter=rp["max_iter"],
                    importance_nested_sampling=rp["do_ins"],
                    outputfiles_basename=rp["outputfiles_basename"],\
                    init_MPI=False)

    #progress.stop()

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    # Run Multinest
    ret = main()
    sys.exit(ret)






