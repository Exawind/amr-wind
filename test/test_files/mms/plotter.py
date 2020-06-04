#!/usr/bin/env python3

import argparse
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

plt.rc("text", usetex=True)
cmap_med = [
    "#F15A60",
    "#7AC36A",
    "#5A9BD4",
    "#FAA75B",
    "#9E67AB",
    "#CE7058",
    "#D77FB4",
    "#737373",
]
cmap = [
    "#EE2E2F",
    "#008C48",
    "#185AA9",
    "#F47D23",
    "#662C91",
    "#A21D21",
    "#B43894",
    "#010202",
]
dashseq = [
    (None, None),
    [10, 5],
    [10, 4, 3, 4],
    [3, 3],
    [10, 4, 3, 4, 3, 4],
    [3, 3],
    [3, 3],
]
markertype = ["s", "d", "o", "p", "h"]


if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(description="A simple plot tool")
    parser.add_argument(
        "-f", "--fdirs", help="Folders to plot", type=str, required=True, nargs="+"
    )
    args = parser.parse_args()

    fields = ["u", "v", "w"]

    lst = []
    for k, fdir in enumerate(args.fdirs):
        df = pd.read_csv(os.path.join(fdir, "mms.log"), delim_whitespace=True)
        df["res"] = float(fdir)
        lst.append(df.iloc[-1])

        for field in fields:
            plt.figure(field)
            p = plt.plot(df.time, df[f"L2_{field}"], lw=2, color=cmap[k], label=fdir)
            p[0].set_dashes(dashseq[k])

    df = pd.DataFrame(lst)
    idx = 1
    theory_order = 2
    df["theory"] = df["L2_u"].iloc[idx] * (df.res.iloc[idx] / df.res) ** theory_order

    fname = "plots.pdf"
    with PdfPages(fname) as pdf:
        for k, field in enumerate(fields):
            plt.figure(field)
            ax = plt.gca()
            plt.xlabel(r"$t$", fontsize=22, fontweight="bold")
            plt.ylabel(f"$L_2({field})$", fontsize=22, fontweight="bold")
            plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
            plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
            legend = ax.legend(loc="best")
            plt.tight_layout()
            pdf.savefig(dpi=300)

            plt.figure(f"ooa")
            plt.loglog(
                df.res,
                df[f"L2_{field}"],
                lw=2,
                color=cmap[k],
                marker=markertype[k],
                mec=cmap[k],
                mfc=cmap[k],
                ms=10,
                label=field,
            )

        plt.figure(f"ooa")
        plt.loglog(df.res, df.theory, lw=1, color=cmap[-1], label="2nd order")
        ax = plt.gca()
        plt.xlabel(r"$N$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$L_2$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)
