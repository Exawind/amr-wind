"""Example script for reading and writing a AMR-Wind plot file."""

import argparse
import pathlib
import amrex.space3d as amr
from amrex_plotfile import AmrexPlotFile
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="An example tool to manipulate plot file data"
    )
    parser.add_argument(
        "-f",
        "--fname",
        help="AMR-Wind plot file",
        required=True,
        type=str,
    )
    args = parser.parse_args()

    amr.initialize([])

    hname = pathlib.Path(args.fname) / "Header"
    wname = pathlib.Path(hname.parent.name)
    plt = AmrexPlotFile(hname)
    mfs = plt()
    plt.write(wname)


if __name__ == "__main__":
    main()
