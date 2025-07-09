# -*- coding: utf-8 -*-

"""\
Utilities from amrex operations
-------------------------------

"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator


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


def ncomp_from_field(field):
    return 3 if field == "velocity" else 1


def interpolate(plt, plti):
    """Interpolate from one plt file to another"""
    assert plt.nlevels == plti.nlevels
    for ilev in range(plt.nlevels):
        assert plt.ngrids[ilev] == plti.ngrids[ilev]
        for igrid in range(plt.ngrids[ilev]):
            original = plt.mfs[ilev].to_xp()[igrid]
            xg, yg, zg = plt.coordinates(ilev, igrid)

            datai = plti.mfs[ilev].to_xp()[igrid]
            xgi, ygi, zgi = plti.coordinates(ilev, igrid)

            for nc in range(plt.ncomp):
                data = original[:, :, :, nc]
                interp = RegularGridInterpolator(
                    (xg[:, 0, 0], yg[0, :, 0], zg[0, 0, :]),
                    data,
                    bounds_error=False,
                    fill_value=None,
                )
                datai = plti.mfs[ilev].to_xp()[0]
                datai[:, :, :, nc] = interp((xgi, ygi, zgi))
