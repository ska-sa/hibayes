#!/usr/bin/env python

"""
# hi_plot.py

Jonathan Zwart
Danny Price
Gianni Bernardi
April 2016

Usage:

./hi_plot.py config_file.ini

"""

import sys
import time
import os
import pprint
from hibayes.parse_config import parse_config
import pylab

from hibayes import contour_plot
from hibayes.utils import *

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print 'usage:'
        print './hi_plot.py config_filename.ini'
        print
        sys.exit(0)

param_file = sys.argv[-1]
# load runtime parameters
rp = parse_config(param_file)


#-------------------------------------------------------------------------------


@profile
def main():
    """
    """

    print 'Settings file is %s' % param_file

    # Insert custom tweaks here e.g.
    #rp["plotRanges"]['C']=[0,200]

    triangle = 'triangle.pdf'
    autoscale = False

    chain = pylab.loadtxt('%s/1-post_equal_weights.dat' % rp["outdir"])

    # Make a budget plot
    bundle = contour_plot.contourTri(chain,
                                     line=True, outfile='%s/%s' % (rp["outdir"], triangle),
                                     col=('red', 'blue'), labels=rp["parameters"],
                                     ranges=rp["plotRanges"],autoscale=autoscale)

    # Plot for publication
    line = False
    autoscale = True
    title = ''
    truth = None
    extn = 'pdf'
    binsize = 25
    bundle = contour_plot.contourTri(chain,
                                     line=line,
                                     outfile='%s/triangle-publn.%s' % (rp["outdir"], extn),
                                     labels=rp["parameters"],
                                     ranges=rp["plotRanges"], truth=truth,
                                     autoscale=autoscale, title=title,
                                     binsize=binsize, labelDict=rp["labelDict"])

    stats = fetchStats(rp["outdir"], rp["parameters"], rp["plotTruth"])
    printLaTeX(rp["parameters"], stats, dump=rp["outdir"])

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    ret = main()
    sys.exit(ret)
