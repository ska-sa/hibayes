"""
# profile_support.py

Tools/decorator for code profiling
"""

import __builtin__

try:
    profile = __builtin__.profile
except AttributeError:
    # No line profiler, provide a pass-through version
    def profile(func):
        return func
