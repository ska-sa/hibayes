"""
# constants.py

Constants used for calculations within ledabayes package
"""

from numpy import pi, sqrt, log

# Constants
muJy2Jy = 1.0e-6
sqDeg2sr = 4.0 * pi * pi / 129600.0
sqrtTwo = sqrt(2.0)
Jy2muJy = 1.0e6
beamFac = pi / (4.0 * log(2.0))