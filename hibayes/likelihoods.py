"""
# likelihoods.py

Place likelihood and prior functions here
"""

import sys
import numpy
from hibayes.priors import Priors
from hibayes.profile_support import profile
from hibayes.parse_config import parse_config
from hibayes.sky_model import generate_simulated_data
from hibayes.spectral_models import T_HI,T_fg,sigma
#from hibayes.sky_model import nu_1
param_file = sys.argv[-1]

rp = parse_config(param_file)

#-------------------------------------------------------------------------------

if True:
    if rp["ledaSpec"] is None:
        Tmeas, freqs = generate_simulated_data(rp)
    else:
        freqs = numpy.genfromtxt(rp["ledaFreqs"])
        Tmeas = numpy.genfromtxt(rp["ledaSpec"])
        FREQ_MIN = freqs[0]
        FREQ_MAX = freqs[-1]

#-------------------------------------------------------------------------------

pri = Priors()
@profile
def logprior(cube, ndim, nparams, fg_only=False, bg_only=False, log_prior=False):
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
        # Comment this out if you're fitting for the HI too
        cube[1]=pri.UniformPrior(cube[1], rp["FREQ_MIN"],rp["FREQ_MAX"]) # nu_HI/MHz
    elif fg_only and not bg_only:
        # Comment this out if you're fitting for the foregrounds only
        cube[1]=pri.GeneralPrior(cube[1],'U', rp["FREQ_MIN"], rp["FREQ_MAX"])

    cube[2] = pri.GeneralPrior(cube[2], 'U', rp["SIGMA_HI_MIN"], rp["SIGMA_HI_MAX"])  # sigma_HI/MHz

    for ic in range(ndim - 3):
        cube[ic + 3] = pri.GeneralPrior(cube[ic + 3], 'U', -rp["BP_PRIOR_RANGE"], +rp["BP_PRIOR_RANGE"])
    
    #print "%2.2f | %2.2f | %2.2f" % (cube[0], cube[1], cube[2])
    return


#-------------------------------------------------------------------------------

@profile
def loglike(cube, ndim, nparams):
    """
    Form of lhood assumes uncertainties are gaussian and uncorrelated
    The latter can be extended later
    
    """
    #numpy.savetxt('Tmeas.txt',Tmeas)	# Test whether reads correct spectrum when run on data
    #numpy.savetxt('freqs.txt',freqs)	# Test whether reads correct frequency when run on data

    A_HI = cube[0]
    nu_HI = cube[1]
    sigma_HI = cube[2]
    #c=[cube[ic+3] for ic in range(ndim)]

    c_fit = [cube[i + 3] for i in range(ndim - 3)]

    loglike = 0.0
    for idatum in range(len(freqs)):
        Tsky = 1.0e-3 * T_HI(A_HI, nu_HI, sigma_HI, freqs[idatum]) + \
          T_fg(rp["nu_1"], c_fit, rp["nc_fit"], freqs[idatum])
        if rp["spectrum_errors"] is not None:
            sig = errors[idatum]
        else:
            sig = sigma(Tsky, rp["BW"], rp["tObs"])
        chisq = 0.5 * ((Tmeas[idatum] - Tsky) / sig) ** 2.0
        prefactor = 0.5 * numpy.log(2.0 * numpy.pi * sig ** 2.0)
        loglike -= prefactor + chisq

    return loglike

#-------------------------------------------------------------------------------
