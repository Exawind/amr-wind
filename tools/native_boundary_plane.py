import argparse
import glob
import pathlib
import amrex.space3d as amr
from amrex_plotfile import AmrexPlotFile
import numpy as np


def parse_fields_from_header_name(hname):
    return hname.split("_")[-1]


def parse_orientations_from_header_name(hname):
    return int(hname.split("_")[-2])


def main():
    parser = argparse.ArgumentParser(description="A tool to refine boundary planes")
    parser.add_argument(
        "-f",
        "--fdir",
        help="AMR-Wind directory with boundary data",
        required=True,
        type=str,
    )
    args = parser.parse_args()

    amr.initialize([])

    for fname in sorted(glob.glob(f"{args.fdir}/bndry_output" + "*")):
        print(f"Reading {fname}")
        headers = [pathlib.Path(x) for x in glob.glob(f"{fname}/Header_*")]
        fpath = pathlib.Path(fname)
        wname = pathlib.Path(fpath.parent.name) / pathlib.Path(fpath.name)

        for hname in headers:
            field = parse_fields_from_header_name(hname.name)
            ori = parse_orientations_from_header_name(hname.name)
            plt = AmrexPlotFile(hname)
            mfs = plt()
            for mf in mfs:
                assert mf.box_array().size == 1
            plt.write(wname)


if __name__ == "__main__":
    main()

# # Read plt file
# #dname = "/Users/mhenryde/exawind/source/amr-wind/build/test/test_files/abl_bndry_output_native"
# dname = "/Users/mhenryde/exawind/source/amr-wind/build/test/test_files/abl_bndry_output_amr_native"
# fname = dname + "/bndry_files/bndry_output00004/Level_0/velocity_0"
# # fname = dname + "/plt00000/Level_0/Cell"
# mf = amr.VisMF.Read(fname)
# mf_xp = mf.to_xp()
# # mf_xp[0] = np.ones(mf_xp[0].shape)
# # fname += '_new'
# # amr.VisMF.Write(mf, fname)
# print(mf_xp[0].shape)
# print(len(mf_xp))
# # Just checking stuff "works"
# small_end = amr.IntVect()
# big_end = amr.IntVect(2, 3, 4)
# b = amr.Box(small_end, big_end)
# print(b)
