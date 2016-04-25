"""
# utils.py

Useful functions
"""

import os

import numpy
import scipy
from scipy import stats
from scipy.interpolate import interp1d
from scipy.special import erf
import pymultinest

from hibayes.profile_support import profile

#-------------------------------------------------------------------------------

# Constants
from .constants import *

#-------------------------------------------------------------------------------

def gaussian(x, mu, sig, norm=True):
    gauss = numpy.exp(-0.5 * ((x - mu) / sig) ** 2)
    if norm:
        gauss *= 1.0 / (sqrtTwo * sqrt(pi) * sig)
    return gauss


#-------------------------------------------------------------------------------

@profile
def dump_variable_values(module, moduleFile, verbose=False):
    """
    http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order
    """

    # Retain only last (active) occurrence of any repeated variable
    mylist = [x for i, x in enumerate(dir(module)) if x not in dir(module)[i + 1:]]

    f = open(moduleFile, 'w')

    for i in range(len(dir(module))):
        line = ['%s %s' % (x, getattr(module, x)) for x in dir(module)][i]
        if verbose: print line
        f.write('%s\n' % line)

    f.close()

    return


#-------------------------------------------------------------------------------

@profile
def peak_confidence(vector, bins=None):
    """
    For a vector, bin data and find peak position
    See http://stackoverflow.com/questions/18818905/find-the-x-value-corresponding-to-a-histogram-max
    """

    if bins is None: bins = 100

    n, b = numpy.histogram(vector, bins=bins)
    bin_max = numpy.where(n == n.max())
    bin_max_posn_lower = bin_max[0][0]
    peak = (b[bin_max_posn_lower] + b[bin_max_posn_lower + 1]) / 2.0

    return peak


#-------------------------------------------------------------------------------

@profile
def calculate_confidence(vector, alpha=0.68, ret_all=False):
    """
    from stacker.py (modified)
    """

    percentile_median = 50.0
    percentile_low = percentile_median - (100.0 * alpha / 2.0)
    percentile_high = percentile_median + (100.0 * alpha / 2.0)

    median = stats.scoreatpercentile(vector, percentile_median)
    err_low = median - stats.scoreatpercentile(vector, percentile_low)
    err_high = stats.scoreatpercentile(vector, percentile_high) - median

    if ret_all:
        return median, err_low, err_high, \
               stats.scoreatpercentile(vector, percentile_low), \
               stats.scoreatpercentile(vector, percentile_high)
    else:
        return median, err_low, err_high


#-------------------------------------------------------------------------------

@profile
def mean_confidence_interval(data, confidence=0.95):
    """
    from http://stackoverflow.com/questions/15033511/compute-a-confidence-interval-from-sample-data
    """
    a = 1.0 * numpy.array(data)
    n = len(a)
    m, se = numpy.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t._ppf((1 + confidence) / 2., n - 1)
    return m, m - h, m + h


#-------------------------------------------------------------------------------

@profile
def medianArray(bins):
    """
    Could use map/reduce/ufuncs etc. to compute the median bins
    But for now just do this manually...
    """
    import numpy
    nbins = len(bins)

    bin_medians = numpy.zeros(nbins - 1)
    for ibin in xrange(nbins - 1):
        posts = (bins[ibin], bins[ibin + 1])
        bin_medians[ibin] = numpy.median(posts)
        #print ibin,bins[ibin],bin_medians[ibin]

    return bin_medians


#-------------------------------------------------------------------------------

@profile
def fetchStats(outdir, parameters, truth):
    """
    """

    print '-> analysing summary stats:'

    n_params = len(parameters)
    statsf = '%s/1-' % outdir

    x = pymultinest.analyse.Analyzer(n_params, outputfiles_basename=statsf)
    y = x.get_stats()
    bf = x.get_best_fit()

    Zg = y['global evidence']
    stats = y['marginals']

    B = bf['log_likelihood']

    summary = {}
    for ip, param in enumerate(parameters):
        s = stats[ip]
        b = bf['parameters'][ip]
        print s['1sigma']
        summary[param] = (b, s['median'], s['1sigma'][0], s['1sigma'][-1])

    # ugliest syntax ever!
    print '\n# truth param bestfit median lower upper'
    for param in parameters:
        print '%7s' % param, '%.2f' % truth[param], ' '.join(['%.2f' % s for s in summary[param]])

    print '****Global log-evidence is %f' % Zg

    return summary


#-------------------------------------------------------------------------------

@profile
def printLaTeX(parameters, statsDict, dump=None):
    """
    ' \\\\\n'.join([' & '.join(map(str,line)) for line in a])
    """

    if dump is not None:
        outf = '%s/params.tex' % dump
        out = open(outf, 'w')

    for ip, param in enumerate(parameters):
        val = statsDict[param]
        line = '%6s & $%5.2f_{%6.3f}^{%6.3f}$ \\\\' % (param, val[0], val[2], val[3])
        if dump is not None:
            out.write('%s\n' % line)
        else:
            print line

    if dump is not None:
        out.close()
        print '\n-> writing summary stats to \input{%s}' % outf

    return

#-------------------------------------------------------------------------------

@profile
def fit_polynomial(x, y, order):
    z = numpy.polyfit(x, y, order)
    p = numpy.poly1d(z)

    return z, p, p(x)

#-------------------------------------------------------------------------------

@profile
def single_gaussian(A, x0, sig, freqs, do_erf=False):
    """
    """
    #f(x)=a*exp(-(x-x0)**2/sigma0**2)+b*exp(-(x-x1)**2/sigma1**2)+c*exp(-(x-x2)**2/sigma2**2)
    model = A * numpy.exp(-0.5 * ((freqs - x0) / sig) ** 2)
    # Prefactor of norm doesn't seem to matter?
    if do_erf:
        norm = sqrt(2.0 * pi * sig ** 2) * \
               0.5 * (erf((x0 - freqs[0]) / (sqrtTwo * sig)) + erf((freqs[-1] - x0) / (sqrtTwo * sig)))
        model *= 1.0 / norm
        #model *= sqrt(2.0*pi*sigma**2)
    #if numpy.max(model) != 0.0: model *= 1.0/numpy.max(model)

    return model

#-------------------------------------------------------------------------------

@profile
def read_bandpass(f='bandpass.txt'):
    """
    """

    bpData = numpy.loadtxt(f)

    freqs = bpData[:, 0]
    bandpass = bpData[:, 1]
    #nfreqs=len(freqs)

    return bandpass, freqs

#-------------------------------------------------------------------------------
