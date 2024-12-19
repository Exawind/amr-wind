"""\
A tool to generate native boundary planes header files
------------------------------------------------------

This script is only for generating the header files for boundary plane
data that does not already have them. Most likely they were generated
with AMR-Wind version 3.2.0 and prior.

"""

import argparse
import glob
import pathlib
import itertools
import amrex.space3d as amr
import amrex_utils as au
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="A tool to generate native boundary planes header files"
    )
    parser.add_argument(
        "-f",
        "--fdir",
        help="AMR-Wind directory with boundary data",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-i",
        "--iname",
        help="AMR-Wind input file that generated the boundary data",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Force overwrite of existing header files",
    )
    args = parser.parse_args()

    amr.initialize([])
    spacedim = amr.Config.spacedim

    spath = pathlib.Path(args.fdir)
    tname = "time.dat"
    time_file = spath / tname
    times = pd.read_csv(time_file, sep="\\s+", names=["step", "time"], header=None)
    pfx = "bndry_output"
    lvl_pfx = "Level_"

    coord_system = 0
    bwidth = 0

    pp = amr.ParmParse("")
    pp.addfile(args.iname)
    pp_geom = amr.ParmParse("geometry")
    prob_lo = [
        pp_geom.get_real("prob_lo", 0),
        pp_geom.get_real("prob_lo", 1),
        pp_geom.get_real("prob_lo", 2),
    ]
    prob_hi = [
        pp_geom.get_real("prob_hi", 0),
        pp_geom.get_real("prob_hi", 1),
        pp_geom.get_real("prob_hi", 2),
    ]
    pp_amr = amr.ParmParse("amr")
    n_cell = [
        pp_amr.get_int("n_cell", 0),
        pp_amr.get_int("n_cell", 1),
        pp_amr.get_int("n_cell", 2),
    ]

    ref_ratio = []
    with open(args.iname, "r") as f:
        for line in f:
            if "amr.ref_ratio" in line:
                line = line.split()
                ref_ratio = [int(x) for x in line[2:]]

    for fname in sorted(glob.glob(f"{args.fdir}/{pfx}" + "*")):
        print(f"Generating Header files for data in {fname}")
        fpath = pathlib.Path(fname)
        step = int(fpath.name.replace(pfx, ""))
        time = (times.time[times.step == step]).values[0]

        finest_level = len(sorted(glob.glob(f"{fname}/{lvl_pfx}" + "*"))) - 1
        nlevels = finest_level + 1
        mf_h_names = sorted(glob.glob(f"{fname}/{lvl_pfx}0/" + "*" + "_H"))

        fields = []
        oris = []
        for mf_h_name in mf_h_names:
            fields.append(pathlib.Path(mf_h_name).name.split("_")[0])
            oris.append(int(pathlib.Path(mf_h_name).name.split("_")[1]))
        fields = list(set(fields))
        oris = list(set(oris))

        cell_sizes = [
            [(prob_hi[i] - prob_lo[i]) / (n_cell[i]) for i in range(spacedim)]
        ]
        for ilev in range(1, nlevels):
            cell_sizes.append(
                [x / (ref_ratio[ilev - 1] ** ilev) for x in cell_sizes[0]]
            )

        level_steps = [step] * nlevels

        for field, ori in itertools.product(fields, oris):
            hname = fpath / f"Header_{ori}_{field}"
            if hname.exists() and not args.overwrite:
                print(f"{hname} exists already. Skipping.")
                continue

            mfs = []
            for ilev in range(nlevels):
                mf_h_name = fpath / f"{lvl_pfx}{ilev}" / f"{field}_{ori}_H"
                mfs.append(amr.VisMF.Read(str(mf_h_name).replace("_H", "")))
            ncomp = mfs[0].num_comp
            vnames = au.variable_names(field, ncomp)

            bas = [mf.box_array() for mf in mfs]

            ngrids = [ba.size for ba in bas]
            for ilev in range(nlevels):
                assert ngrids[ilev] == 1

            normal = au.normal_from_ori(ori)
            glohis = [[] for _ in range(nlevels)]
            for ilev in range(nlevels):
                nglohis = []
                for igrid in range(ngrids[ilev]):
                    glos = []
                    ghis = []
                    for idim in range(spacedim):
                        if au.is_low(ori):
                            glo = (
                                prob_lo[idim]
                                if idim != normal
                                else prob_lo[idim] - cell_sizes[ilev][idim]
                            )
                            ghi = (
                                prob_hi[idim]
                                if idim != normal
                                else prob_lo[idim] + cell_sizes[ilev][idim]
                            )
                        else:
                            glo = (
                                prob_lo[idim]
                                if idim != normal
                                else prob_hi[idim] - cell_sizes[ilev][idim]
                            )
                            ghi = (
                                prob_hi[idim]
                                if idim != normal
                                else prob_hi[idim] + cell_sizes[ilev][idim]
                            )
                        glos.append(glo)
                        ghis.append(ghi)
                    nglohis.append([glos, ghis])
                glohis[ilev] = nglohis

            with open(hname, "w") as f:
                f.write("HyperCLaw-V1.1\n")
                f.write(f"{ncomp}\n")
                f.write("\n".join(vnames))
                f.write("\n")
                f.write(f"{spacedim}\n")
                f.write(f"{time:.17g}\n")
                f.write(f"{finest_level}\n")
                f.write(" ".join([f"{x:.17g}" for x in prob_lo]))
                f.write(" \n")
                f.write(" ".join([f"{x:.17g}" for x in prob_hi]))
                f.write(" \n")
                f.write(" ".join([f"{x:.17g}" for x in ref_ratio]))
                f.write(" \n")

                line = "("
                for ilev in range(nlevels):
                    ba = bas[ilev]
                    small_end = ba[0].small_end
                    big_end = ba[0].big_end
                    line += "("
                    for idim in range(spacedim):
                        end = "," if idim != spacedim - 1 else ")"
                        line += str(small_end[idim]) + end
                    line += " ("
                    for idim in range(spacedim):
                        end = "," if idim != spacedim - 1 else ")"
                        line += str(big_end[idim]) + end
                    line += " ("
                    for idim in range(spacedim):
                        end = "," if idim != spacedim - 1 else ")"
                        line += str(0) + end
                    end = " (" if ilev != finest_level else ""
                    line += ")" + end
                f.write(f"{line} \n")

                f.write(" ".join([str(x) for x in level_steps]))
                f.write(" \n")

                for ilev in range(nlevels):
                    f.write(" ".join([str(x) for x in cell_sizes[ilev]]))
                    f.write(" \n")

                f.write(f"{coord_system}\n")
                f.write(f"{bwidth}\n")

                mf_names = [f"{lvl_pfx}{x}/{field}_{ori}" for x in range(nlevels)]
                for ilev in range(nlevels):
                    f.write(f"{ilev} {ngrids[ilev]} {time:.17g}\n")
                    f.write(f"{level_steps[ilev]}\n")
                    for igrid in range(ngrids[ilev]):
                        for idim in range(spacedim):
                            glo = glohis[ilev][igrid][0][idim]
                            ghi = glohis[ilev][igrid][1][idim]
                            f.write(f"{glo:.17g} {ghi:.17g}\n")
                    f.write(f"{mf_names[ilev]}\n")


if __name__ == "__main__":
    main()
