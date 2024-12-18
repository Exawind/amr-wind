# -*- coding: utf-8 -*-

"""\
Utilities from amrex operations
-------------------------------

"""

import numpy as np


def is_low(ori):
    """Return whether the orientation is a low face."""
    if ori == 0 or ori == 1 or ori == 2:
        return True
    else:
        return False


def normal_from_ori(ori):
    """Return the normal directions given an orientation."""
    if ori == 0 or ori == 3:
        return 0
    elif ori == 1 or ori == 4:
        return 1
    elif ori == 2 or ori == 5:
        return 2
    else:
        raise Exception("Invalid ori")
    return -1


def perpendicular_from_ori(ori):
    """Return the perpendicular directions given an orientation."""
    normal = normal_from_ori(ori)
    if normal == 0:
        return [1, 2]
    elif normal == 1:
        return [0, 2]
    elif normal == 2:
        return [0, 1]
    else:
        raise Exception("Invalid normal")
    return [-1, -1]


def variable_names(field, ncomp):
    if ncomp == 1:
        return [field]
    elif ncomp == 3:
        return [f"{field}{x}" for x in ["x", "y", "z"]]
    else:
        raise Exception("Invalid ncomp")


def slice_from_normal(normal, i, comp):
    """Returns a numpy slice for convenient indexing in boundary planes"""
    if normal == 0:
        return np.s_[i, :, :, comp]
    if normal == 1:
        return np.s_[:, i, :, comp]
    if normal == 2:
        return np.s_[:, :, i, comp]
    else:
        raise Exception("Invalid normal")
