#!/usr/bin/env python3

"""
Python interface for libdom.
"""

# libdom interface
import numpy as np
from numpy.ctypeslib import ndpointer
import ctypes
import os

_doublepp = ndpointer(dtype=np.uintp, ndim=1, flags='C')

_dll = ctypes.CDLL(os.path.join(os.environ['MODCELLHPC_PATH'],'io','libdom','libdom.so'))

_fdom = _dll.fdom
_fdom.argtypes = [ctypes.c_int, ctypes.c_int, _doublepp, ctypes.c_void_p]
_fdom.restype = None

def fdom(x):
    """ Find non-dominated rows

    Parameters
    ----------
    x : m by n array of doubles
        Represents Pareto front.

    Returns
    -------
     y : m by 1 array
        Contains True if a row is non-dominated and False otherwise
    """

    # Make sure x is double and contigous:
    x = np.ascontiguousarray(x, dtype=np.double)
    # Create y
    y = np.zeros((x.shape[0],1), dtype=x.dtype) # Should be int
    # run fdom with appropriate parameters
    xpp = (x.__array_interface__['data'][0] + np.arange(x.shape[0])*x.strides[0]).astype(np.uintp)
    yp = ctypes.c_void_p(y.ctypes.data)
    m = ctypes.c_int(x.shape[0])
    n = ctypes.c_int(x.shape[1])
    _fdom(m, n, xpp, yp)
    # convert y to bool
    y = y.astype(dtype=bool)
    return np.logical_not(y) # y is negated since the C implementation returns True if a solution is dominated


