# coding: utf-8
import argparse

import h5py
import hdf5plugin
import numpy as np


def modify_attributes(fname):
    """Modifies HDF5 attributes from amrex plt files to be compatible with yt"""
    hf = h5py.File(fname, "r+")
    attr_keys = ["dim", "num_levels", "finest_level"]
    for key in attr_keys:
        if np.ndim(hf.attrs[key]):
            hf.attrs.__setitem__(key, hf.attrs[key].item())
    chombo_global_keys = ["SpaceDim"]
    for key in chombo_global_keys:
        if np.ndim(hf["Chombo_global"].attrs[key]):
            hf["Chombo_global"].attrs.__setitem__(
                key, hf["Chombo_global"].attrs[key].item()
            )
    hf.close()


def main():
    """Plot data."""
    parser = argparse.ArgumentParser(
        description="Modifies HDF5 attributes from amrex plt files to be compatible with yt"
    )
    parser.add_argument(
        "-f",
        "--fname",
        help="plt h5 files to replace attrs in",
        required=True,
        type=str,
        nargs="+",
    )
    args = parser.parse_args()

    for fname in args.fname:
        modify_attributes(fname)


if __name__ == "__main__":
    main()
