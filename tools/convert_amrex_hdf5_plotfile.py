import h5py
import hdf5plugin
import numpy as np
from pathlib import Path
import argparse


def main():
    """Convert an AMReX HDF5 plotfile 3D npy file."""
    parser = argparse.ArgumentParser(
        description="Convert an AMReX HDF5 plotfile 3D npy file"
    )
    parser.add_argument(
        "-f",
        "--fname",
        help="AMReX HDF5 plotfile name",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-o",
        "--odir",
        help="Directory for output files",
        default=".",
        type=str,
    )
    args = parser.parse_args()

    with h5py.File(args.fname, mode="r") as f:
        dct = {}
        for key in f.attrs.keys():
            if "component_" in key:
                field = f.attrs[key].decode("UTF-8")
                component = int(key.split("_")[1])
                dct[field] = component
        nvars = len(dct)

        # Characterize box array
        ba_shape0 = f["level_0"]["boxes"][0]
        ba_di = ba_shape0[3] - ba_shape0[0] + 1
        ba_dj = ba_shape0[4] - ba_shape0[1] + 1
        ba_dk = ba_shape0[5] - ba_shape0[2] + 1
        ba_var_len = ba_di * ba_dj * ba_dk
        boxstride = ba_var_len * nvars
        nba = len(f["level_0"]["boxes"])

        # Check that di/dj/dk are true for all shapes
        for box in f["level_0"]["boxes"]:
            assert box[0] + ba_di - 1 == box[3], box
            assert box[1] + ba_dj - 1 == box[4], box
            assert box[2] + ba_dk - 1 == box[5], box

        # Characterize domain
        dom_nx = f["level_0"].attrs["prob_domain"][3] + 1
        dom_ny = f["level_0"].attrs["prob_domain"][4] + 1
        dom_nz = f["level_0"].attrs["prob_domain"][5] + 1

        # Reformat 1D into 3D data, iterating over variables
        data = np.zeros((nvars, dom_nx, dom_ny, dom_nz))
        for i, (field, comp) in enumerate(dct.items()):
            var_offset = comp * ba_var_len

            for lb in range(nba):
                ba_inds = f["level_0"]["boxes"][lb]
                lo_i, lo_j, lo_k = ba_inds[0], ba_inds[1], ba_inds[2]
                hi_i, hi_j, hi_k = ba_inds[3] + 1, ba_inds[4] + 1, ba_inds[5] + 1

                hdf5_data = f["level_0"]["data:datatype=0"][
                    var_offset
                    + lb * boxstride : var_offset
                    + lb * boxstride
                    + ba_var_len
                ]

                # Reshape into 3D
                np_data = np.reshape(
                    hdf5_data, (hi_i - lo_i, hi_j - lo_j, hi_k - lo_k), order="F"
                )

                # Place little 3D data into the big domain
                data[i, lo_i:hi_i, lo_j:hi_j, lo_k:hi_k] = np_data

        savename = f"{args.fname}.npy"
        np.save(Path(args.odir, savename), data)


if __name__ == "__main__":
    main()
