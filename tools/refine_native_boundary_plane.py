import argparse
import glob
import pathlib
import shutil
import amrex.space3d as amr
from amrex_plotfile import AmrexPlotFile
import amrex_utils as au
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def parse_fields_from_header_name(hname):
    return hname.split("_")[-1]


def parse_orientations_from_header_name(hname):
    return int(hname.split("_")[-2])


def refine(plt, ori, refinement_ratio):
    """Update data in a boundary plane for integer refinement."""
    normal = au.normal_from_ori(ori)
    perp = au.perpendicular_from_ori(ori)

    for ilev in range(plt.nlevels):
        small_endi = plt.prob_domain[ilev].small_end
        big_endi = plt.prob_domain[ilev].big_end
        for p in perp:
            small_endi[p] *= refinement_ratio
            big_endi[p] = (big_endi[p] + 1) * refinement_ratio - 1
        plt.prob_domain[ilev] = amr.Box(small_endi, big_endi)

        plt.cell_sizes[ilev] = [x / refinement_ratio for x in plt.cell_sizes[ilev]]

        dx = plt.cell_sizes[ilev][normal]
        for igrid in range(plt.ngrids[ilev]):
            glo = plt.glohis[ilev][igrid][0][normal]
            ghi = plt.glohis[ilev][igrid][1][normal]
            offset = (ghi + glo) / 2
            plt.glohis[ilev][igrid][0][normal] = offset - dx
            plt.glohis[ilev][igrid][1][normal] = offset + dx

        bac = plt.mfs[ilev].box_array()
        rr = amr.IntVect(refinement_ratio)
        rr[normal] = 1
        baf = bac.refine(rr)

        dm = plt.mfs[ilev].dm()
        plt.mfs[ilev] = amr.MultiFab(baf, dm, plt.ncomp, 0)
        plt.mfs[ilev].set_val(0.0)


def interpolate(plt, plti, ori, refinement_ratio):
    """Interpolate from one boundary file to another"""
    normal = au.normal_from_ori(ori)
    perp = au.perpendicular_from_ori(ori)

    for ilev in range(plt.nlevels):
        assert plt.mfs[ilev].box_array().size == 1
        assert plt.ngrids[ilev] == 1
        original = plt.mfs[ilev].to_xp()[0]
        los = plt.glohis[ilev][0][0]
        his = plt.glohis[ilev][0][1]
        xlo = los[perp[0]]
        xhi = his[perp[0]]
        ylo = los[perp[1]]
        yhi = his[perp[1]]
        nx = plt.mfs[ilev].box_array()[0].size[perp[0]]
        ny = plt.mfs[ilev].box_array()[0].size[perp[1]]
        dx = plt.cell_sizes[ilev][perp[0]]
        dy = plt.cell_sizes[ilev][perp[1]]
        x = np.linspace(xlo + 0.5 * dx, xhi - 0.5 * dx, nx)
        y = np.linspace(ylo + 0.5 * dy, yhi - 0.5 * dy, ny)

        nxi = plti.mfs[ilev].box_array()[0].size[perp[0]]
        nyi = plti.mfs[ilev].box_array()[0].size[perp[1]]
        dxi = plti.cell_sizes[ilev][perp[0]]
        dyi = plti.cell_sizes[ilev][perp[1]]
        xi = np.linspace(xlo + 0.5 * dxi, xhi - 0.5 * dxi, nxi)
        yi = np.linspace(ylo + 0.5 * dyi, yhi - 0.5 * dyi, nyi)
        xgi, ygi = np.meshgrid(xi, yi, indexing="ij")

        for i in range(original.shape[normal]):
            for nc in range(plt.ncomp):
                slc = au.slice_from_normal(normal, i, nc)
                data = original[slc]
                interp = RegularGridInterpolator(
                    (x, y), data, bounds_error=False, fill_value=None
                )
                datai = plti.mfs[ilev].to_xp()[0]
                datai[slc] = interp((xgi, ygi))


def main():
    parser = argparse.ArgumentParser(description="A tool to refine boundary planes")
    parser.add_argument(
        "-f",
        "--fdir",
        help="AMR-Wind directory with boundary data",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-r",
        "--refinement_ratio",
        help="Refinement ratio",
        default=2,
        type=int,
    )
    args = parser.parse_args()

    amr.initialize([])

    for fname in sorted(glob.glob(f"{args.fdir}/bndry_output" + "*")):
        print(f"Refining {fname}")
        headers = [pathlib.Path(x) for x in glob.glob(f"{fname}/Header_*")]
        fpath = pathlib.Path(fname)
        wname = pathlib.Path(fpath.parent.name) / pathlib.Path(fpath.name)

        for hname in headers:
            field = parse_fields_from_header_name(hname.name)
            ori = parse_orientations_from_header_name(hname.name)
            plt = AmrexPlotFile(hname)
            mfs = plt()

            plti = AmrexPlotFile(hname)
            _ = plti()
            refine(plti, ori, args.refinement_ratio)
            interpolate(plt, plti, ori, args.refinement_ratio)
            plti.write(wname)

    # Copy the time file
    tname = "time.dat"
    spath = pathlib.Path(args.fdir)
    src_time_file = spath / tname
    dst_time_file = pathlib.Path(spath.name) / tname
    shutil.copyfile(src_time_file, dst_time_file)


if __name__ == "__main__":
    main()
