"""
# profile_support.py

Tools/decorator for code profiling
"""

import six
if six.PY2:
    import __builtin__ as builtins
else:
    import builtins

try:
    profile = builtins.profile
except AttributeError:
    # No line profiler, provide a pass-through version
    def profile(func):
        return func
