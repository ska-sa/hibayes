#!/usr/bin/env python

"""
# hi_showfit.py

Reconstruct the (measured) spectrum from Multinest output
Generates an ascii file recon_stats.txt that you can use for plotting

Jonathan Zwart
Danny Price
Gianni Bernardi
April 2016

Usage:

./hi_recon.py config_file.ini

"""

import os
import sys


import numpy as np
import pylab as plt

from hibayes.utils import fetchStats
from hibayes.parse_config import parse_config
from hibayes.spectral_models import T_HI, T_fg

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('usage:')
        print('./hi_recon.py config_filename.ini')
        print()
        sys.exit(0)

param_file = sys.argv[-1]
# load runtime parameters
rp = parse_config(param_file)

def main(c0='#333333', c1='#cc0000'):
    """
    """
    s = fetchStats(rp["outdir"], rp["parameters"], rp["plotTruth"])
    nc = len(s) - 3

    f = np.genfromtxt(rp['ledaFreqs'])
    d = np.genfromtxt(rp['ledaSpec'])
    e = np.genfromtxt(rp['spectrum_errors'])

    T_21f = T_HI(s['A_HI'][0], s['NU_HI'][0], s['SIGMA_HI'][0], f, norm=False, erfs=False) * 1e-3
    fit_coeffs = [s['p%i'%ii][0] for ii in range(len(s) - 3)]
    T_fgf  = T_fg(rp['nu_1'], fit_coeffs, nc, f)

    plt.subplot(3,1,1)
    plt.plot(f, d, label='data', c=c0)
    plt.plot(f, T_fgf+T_21f, label='fit', c=c1)
    plt.ylabel("Temperature [K]")
    plt.legend()

    plt.subplot(3,1,2)
    plt.plot(f, d-(T_fgf+T_21f), label='data - fit', c=c0)
    plt.ylabel("Temperature [K]")
    plt.legend()

    plt.subplot(3,1,3)
    plt.plot(f, d-T_fgf, label='data - foregrounds', c=c0)
    plt.plot(f, T_21f, label='$T_{21}$ fit', c=c1)
    plt.legend()
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Temperature [K]")
    plt.show()


if __name__ == '__main__':
    ret = main()
    sys.exit(ret)
