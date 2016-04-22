#!/usr/bin/env python

"""
# hi_recon.py

Reconstruct the (measured) spectrum from Multinest output

Jonathan Zwart
September 2015

Usage:

./hi_recon.py config_file.ini

"""

import os
import sys


import numpy
from scipy import stats
import pymultinest

from hibayes.utils import calculate_confidence, peak_confidence
from hibayes.spectral_models import T_fg, T_HI
from hibayes.sky_model import nu_1, generate_simulated_data
from hibayes.parse_config import parse_config

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print 'usage:'
        print './hi_plot.py config_filename.ini'
        print
        sys.exit(0)

param_file = sys.argv[-1]
# load runtime parameters
rp = parse_config(param_file)


def func(fr, drawmap=None, fr_1=None, subtractValue=None):
    recon = 1.0e-3 * T_HI(drawmap[0], drawmap[1], drawmap[2], fr) + T_fg(fr_1, drawmap[3:], len(drawmap[3:]), fr)
    if subtractValue is not None:
        recon -= subtractValue
    return recon
    # ff=numpy.vectorize(func,excluded=('drawmap','fr_1'))


# -------------------------------------------------------------------------------

def main():
    """
    """

    # Import the settings variables
    print 'Settings file is %s' % param_file

    # Import the settings variables
    #set_module = importlib.import_module(settingsf)
    #globals().update(set_module.__dict__)

    # Set up the experiment
    #expt=countModel(modelFamily,nlaws,settingsf,dataset,floatNoise)

    f = '%spost_equal_weights.dat' % rp["outstem"]
    f = os.path.join(rp["outdir"], f)

    if rp["ledaSpec"] is None:
        data, freqs = generate_simulated_data(rp)
    else:
        freqs = numpy.genfromtxt(rp["ledaFreqs"])
        data = numpy.genfromtxt(rp["ledaSpec"])

    nfreqs = len(freqs)


    # Load equally-weighted posterior samples
    x = numpy.genfromtxt(f)
    nsamp = x.shape[0]
    ncols = x.shape[1]  # The fifth [seventh] column is the posterior value
    # There must be a better way, but:
    z = numpy.zeros((nsamp, ncols - 1 + nfreqs))
    z[:, :-(nfreqs - 1)] = x
    # Shift posterior values to end
    z[:, -1] = z[:, ncols - 1]  # Copy...
    z[:, ncols - 1] = 0.0  # ...and blank

    # Fetch best-fit parameters and calculate best-fit line
    ana = pymultinest.analyse.Analyzer(ncols - 1, \
                                       outputfiles_basename=os.path.join(rp["outdir"], rp["outstem"]))
    drawml = ana.get_best_fit()['parameters']

    summf = os.path.join(rp["outdir"], '1-summary.txt')
    print summf,ncols
    print numpy.genfromtxt(summf),numpy.genfromtxt(summf).shape
    summary = numpy.genfromtxt(summf)[-1, :]
    drawmap = summary[-(ncols + 1):-2]

    if False:
        print '--> Calculating *ML* reconstruction'
        drawmap = drawml

    # Convert drawmap into correct units etc.
    ymap = numpy.zeros(len(freqs))
    for ifreq, freq in enumerate(freqs):
        ymap[ifreq] = func(freq, drawmap=drawmap, fr_1=nu_1, subtractValue=data[ifreq])
    #ymap=ff(freqs,drawmap=drawmap,fr_1=nu_1)

    for isamp in xrange(nsamp):
        #z[isamp,ncols-1:]=ff(freqs,drawmap=z[isamp,:],fr_1=nu_1)
        for ifreq, freq in enumerate(freqs):
            z[isamp, ncols - 1 + ifreq] = func(freq, drawmap=z[isamp, :], fr_1=nu_1, subtractValue=data[ifreq])

    # Blanking, 0.0 -> NaN
    z[numpy.where(z == 0.0)] = 'NaN'

    # Save the raw reconstructions
    reconf = 'recon_raw.txt'
    reconf = os.path.join(rp["outdir"], reconf)
    recons = z[:, ncols - 1:]
    numpy.savetxt(reconf, recons)

    # Generate stats here...
    s = numpy.zeros((nfreqs, 6))
    s[:, 0] = freqs

    print '# ibin flux fit low high dlower dupper skew kurtosis'
    for ibin in xrange(nfreqs):
        x = recons[:, ibin]
        # Remove NaNs from stats vectors
        # http://stackoverflow.com/questions/11620914/removing-nan-values-from-an-array
        x = x[~numpy.isnan(x)]
        #ss=stats.bayes_mvs(x,alpha=0.68)[0]
        #x*=numpy.power(s[ibin,0]/1.0e6,2.5)

        try:
            ss = numpy.zeros(3)
            ss[0], dlow, dhigh, ss[1], ss[2] = calculate_confidence(x, alpha=0.68, ret_all=True)
        except:
            ss = numpy.nan * numpy.ones(3)
        tt = peak_confidence(x, bins=10)
        #ss*=numpy.power(s[ibin,0]/1.0e6,2.5)
        #print ss[0],tt
        #        s[ibin,1]=ss[0]  # median
        #s[ibin,1]=tt     # peak
        s[ibin, 1] = ymap[ibin]  # MAP
        #print ymap
        s[ibin, 2] = ss[0] - ss[1]  # lower
        s[ibin, 3] = ss[2] - ss[0]  # upper
        s[ibin, 4] = stats.skew(x)  # skewness
        s[ibin, 5] = stats.kurtosis(x)  # kurtosis
        print ibin, s[ibin, 0], s[ibin, 1], dlow, dhigh, ss[1], ss[2], s[ibin, 4], s[ibin, 5]  #,stats.skewtest(x)

    # ...and output to file
    rstatsf = 'recon_stats.txt'
    rstatsf = os.path.join(rp["outdir"], rstatsf)
    hdr = '# freq_MHz T T_lower T_upper skewness kurtosis'
    fid = open(rstatsf, 'w')
    print hdr
    print s

    print '-> Writing reconstructed spectrum to %s' % rstatsf
    fid.write('%s\n' % hdr)
    numpy.savetxt(fid, s)
    fid.close()
    print 'Finished.'

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    ret = main()
    sys.exit(ret)
