#!/usr/bin/env python3

"""
Tests libdom
"""

import numpy as np
from libdomi import fdom

def case(case_n, x, y_exp):
    print (f"Case {case_n}:...")
    print("x:\n",x)
    y = fdom(x)
    print("y:\n",y)
    assert np.alltrue(y == y_exp)

def test():
    x = np.ones((10,5))
    y_exp = np.ones((x.shape[0],1))
    case(1,x, y_exp)

    x = np.array([
        [1, 2.5, 0, 0],
        [3, 0 , 0, 0],
        [1, 2.5, 0, 0],
        [4, 0 , 0, 1]
        ])
    y_exp = np.ones((x.shape[0],1))
    y_exp[1] = False
    case(2,x, y_exp)

    x = np.array([
        [4, 0.5 , 0, 1],
        [1, 2.5, 0, 0],
        [3, 0 , 0, 0],
        [1, 2.5, 0, 0],
        [4, 0 , 0, 1]
        ])
    y_exp = np.ones((x.shape[0],1))
    y_exp[[2,4]] = False

    case(3,x, y_exp)



test()
