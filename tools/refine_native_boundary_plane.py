import argparse
import glob
import pathlib
import shutil
import amrex.space3d as amr
from amrex_plotfile import AmrexPlotFile
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def parse_fields_from_header_name(hname):
    return hname.split("_")[-1]


def parse_orientations_from_header_name(hname):
    return int(hname.split("_")[-2])


def is_low(ori):
    if ori == 0 or ori == 1 or ori == 2:
        return True
    else:
        return False


def normal_from_ori(ori):
    """Return the normal directions given an orientation"""
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
    """Return the perpendicular directions given an orientation"""
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


def refine(plt, ori, refinement_ratio):
    """Update data in a boundary plane for integer refinement."""
    normal = normal_from_ori(ori)
    perp = perpendicular_from_ori(ori)

    for ilev in range(plt.nlevels + 1):
        small_endi = plt.prob_domain[ilev].small_end
        big_endi = plt.prob_domain[ilev].big_end
        for p in perp:
            small_endi[p] *= refinement_ratio
            big_endi[p] = big_endi[p] * refinement_ratio + 1
        plt.prob_domain[ilev] = amr.Box(small_endi, big_endi)

        plt.cell_sizes[ilev] = [x / refinement_ratio for x in plt.cell_sizes[ilev]]

        dx = plt.cell_sizes[ilev][normal]
        for igrid in range(plt.ngrids[ilev]):
            if is_low(ori):
                plo = plt.prob_lo[normal] + refinement_ratio * dx
                plt.glohis[ilev][igrid][0][normal] = plo - dx
                plt.glohis[ilev][igrid][1][normal] = plo + dx
            else:
                phi = plt.prob_hi[normal] - refinement_ratio * dx
                plt.glohis[ilev][igrid][0][normal] = phi - dx
                plt.glohis[ilev][igrid][1][normal] = phi + dx

        ba = amr.BoxArray(plt.prob_domain[ilev])
        dm = plt.mfs[ilev].dm()
        plt.mfs[ilev] = amr.MultiFab(ba, dm, plt.ncomp, 0)
        plt.mfs[ilev].set_val(0.0)


def interpolate(plt, plti, ori, refinement_ratio):
    """Interpolate from one boundary file to another"""
    normal = normal_from_ori(ori)
    perp = perpendicular_from_ori(ori)

    for ilev in range(plt.nlevels + 1):
        assert plt.mfs[ilev].box_array().size == 1
        assert plt.ngrids[ilev] == 1
        original = plt.mfs[ilev].to_xp()[0]
        los = plt.glohis[ilev][0][0]
        his = plt.glohis[ilev][0][1]
        xlo = los[perp[0]]
        xhi = his[perp[0]]
        ylo = los[perp[1]]
        yhi = his[perp[1]]
        nx = plt.prob_domain[ilev].size[perp[0]]
        ny = plt.prob_domain[ilev].size[perp[1]]
        dx = plt.cell_sizes[ilev][perp[0]]
        dy = plt.cell_sizes[ilev][perp[1]]
        x = np.linspace(xlo + 0.5 * dx, xhi - 0.5 * dx, nx)
        y = np.linspace(ylo + 0.5 * dy, yhi - 0.5 * dy, nx)

        nxi = plti.prob_domain[ilev].size[perp[0]]
        nyi = plti.prob_domain[ilev].size[perp[1]]
        dxi = plti.cell_sizes[ilev][perp[0]]
        dyi = plti.cell_sizes[ilev][perp[1]]
        xi = np.linspace(xlo + 0.5 * dxi, xhi - 0.5 * dxi, nxi)
        yi = np.linspace(ylo + 0.5 * dyi, yhi - 0.5 * dyi, nyi)
        xgi, ygi = np.meshgrid(xi, yi, indexing="ij")

        for i in range(original.shape[normal]):
            for nc in range(plt.ncomp):
                slc = slice_from_normal(normal, i, nc)
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
