"""Plot how the RANSConvergence monitor behaves on a KLAxell ABL run.

Reads two things written by the run:

  post_processing/<label>00000.txt   the monitor's own diagnostics
  post_processing/sampling*[0-9].txt raw probe data at the same points

and produces a three-panel figure showing the signal the monitor watches, the
envelope test it applies, and why a difference between successive samples is
not a usable substitute for that test.

Usage: python3 rans_convergence_plot.py <run_directory> [output.png]
"""

import glob
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Timestep of the demonstration run, used to turn a step index into a time
DT = 5.0


def read_monitor(path):
    """Read the monitor's diagnostics file into a dict of columns."""
    data = np.loadtxt(path, skiprows=1)
    # Rows written while the window is still filling are dropped: a window
    # holding a couple of samples has a small spread whatever the flow is
    # doing, so those rows say nothing about convergence
    data = data[data[:, 2] > 0.5]
    return {
        "time": data[:, 0],
        "num_samples": data[:, 1],
        "num_converged": data[:, 3],
        "num_points": data[:, 4],
        "vel_spread": data[:, 6],
        "vel_tol": data[:, 7],
        "tke_spread": data[:, 9],
        "tke_tol": data[:, 10],
        "hold_elapsed": data[:, 11],
        "eta": data[:, 12],
    }


def read_probes(post_dir):
    """Read the ascii probe files into times and per-point arrays.

    Each row of a probe file is x y z, some particle identifiers, and then the
    sampled fields in the order they were requested: u, v, w, tke. Only the
    trailing four columns are relied on, since the identifier block is an
    implementation detail of the particle container.
    """
    files = sorted(glob.glob(os.path.join(post_dir, "sampling*[0-9].txt")))
    times, rows = [], []
    for fname in files:
        step = int(re.search(r"sampling(\d+)\.txt", fname).group(1))
        vals = []
        with open(fname) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 7:
                    continue  # header lines
                vals.append([float(p) for p in parts[-4:]])
        if not vals:
            continue
        times.append(step * DT)
        rows.append(vals)

    order = np.argsort(times)
    times = np.array(times)[order]
    # Shape: (nsamples, npoints, 4) for u, v, w, tke
    return times, np.array(rows)[order]


def fit_decay(times, spreads, threshold=1.0):
    """Least squares of ln(s) against t, the fit the monitor itself performs."""
    mask = spreads > 0.0
    if mask.sum() < 5:
        return None
    slope, intercept = np.polyfit(times[mask], np.log(spreads[mask]), 1)
    rate, amplitude = -slope, np.exp(intercept)
    if rate <= 0.0:
        return None
    resid = np.log(spreads[mask]) - (intercept + slope * times[mask])
    ss_tot = np.sum((np.log(spreads[mask]) - np.mean(np.log(spreads[mask]))) ** 2)
    rsq = 1.0 - (np.sum(resid**2) / ss_tot) if ss_tot > 0 else 0.0
    return {
        "rate": rate,
        "amplitude": amplitude,
        "rsq": rsq,
        "t_cross": np.log(amplitude / threshold) / rate,
    }


def main():
    run_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    out_png = sys.argv[2] if len(sys.argv) > 2 else "rans_convergence_monitor.png"
    post_dir = os.path.join(run_dir, "post_processing")

    monitor_files = glob.glob(os.path.join(post_dir, "convergence*.txt"))
    if not monitor_files:
        sys.exit("No convergence*.txt found in " + post_dir)
    mon = read_monitor(sorted(monitor_files)[0])
    ptimes, pvals = read_probes(post_dir)

    hours = 3600.0
    speed = np.hypot(pvals[:, :, 0], pvals[:, :, 1])
    npts = speed.shape[1]

    # What the criterion actually tests: spread divided by that point's own
    # tolerance, so convergence is the single value one for both quantities
    vel_ratio = mon["vel_spread"] / mon["vel_tol"]
    tke_ratio = mon["tke_spread"] / mon["tke_tol"]

    fig, axes = plt.subplots(3, 1, figsize=(9.5, 10.5), sharex=True)

    # Panel 1: the signal the monitor watches
    ax = axes[0]
    for i in range(npts):
        ax.plot(ptimes / hours, speed[:, i], lw=1.2, label=f"point {i}")
    ax.set_ylabel("horizontal speed  $\\sqrt{u^2+v^2}$  [m/s]")
    ax.set_title(
        "KLAxell ABL approaching pseudo-steady state through a damped "
        "inertial oscillation",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="lower right", ncol=npts)
    ax.grid(alpha=0.3)

    # Panel 2: the envelope test, and the extrapolation
    ax = axes[1]
    ax.semilogy(mon["time"] / hours, vel_ratio, lw=1.4, label="speed envelope")
    ax.semilogy(mon["time"] / hours, tke_ratio, lw=1.4, label="tke envelope")
    ax.axhline(
        1.0, color="k", ls="--", lw=1.2, label="tolerance (converged below)"
    )

    fit = fit_decay(mon["time"], vel_ratio)
    if fit is not None:
        tspan = np.linspace(mon["time"][0], max(fit["t_cross"], mon["time"][-1]), 200)
        ax.semilogy(
            tspan / hours,
            fit["amplitude"] * np.exp(-fit["rate"] * tspan),
            ":",
            color="C0",
            lw=1.6,
            label=(
                f"exponential fit, $R^2$={fit['rsq']:.2f}, "
                f"crossing at {fit['t_cross'] / hours:.1f} h"
            ),
        )
    # Mark the stretches where every point was within tolerance but the
    # criterion had not yet held long enough to stop the run
    holding = mon["hold_elapsed"] > 0.0
    if holding.any():
        ax.fill_between(
            mon["time"] / hours,
            *ax.get_ylim(),
            where=holding,
            color="C2",
            alpha=0.18,
            label="all points within tolerance, holding",
        )

        # Point at the dip, which is the whole reason the hold exists
        dip = mon["time"][np.argmax(holding)] / hours
        ax.annotate(
            "envelope dips below tolerance\nat a turning point;\n"
            "rejected by the hold requirement",
            xy=(dip, 1.0),
            xytext=(dip - 6.0, 0.42),
            fontsize=8.5,
            arrowprops={"arrowstyle": "->", "lw": 1.0, "color": "C2"},
        )

    ax.set_ylabel("spread / tolerance")
    ax.set_title(
        "What the monitor tests: peak-to-trough spread over a trailing window",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3, which="both")

    # Panel 3: why a successive-sample difference is not a substitute
    ax = axes[2]
    worst = int(np.argmax(np.ptp(speed, axis=0)))
    naive = np.abs(np.diff(speed[:, worst]))
    ax.semilogy(
        ptimes[1:] / hours,
        np.maximum(naive, 1e-12),
        lw=1.0,
        color="C3",
        label=f"|change between samples|, point {worst}",
    )
    ax.semilogy(
        mon["time"] / hours,
        mon["vel_spread"],
        lw=1.6,
        color="C0",
        label="window spread, worst point",
    )
    ax.axhline(
        mon["vel_tol"][0],
        color="k",
        ls="--",
        lw=1.2,
        label="speed tolerance",
    )
    ax.set_ylabel("speed change [m/s]")
    ax.set_xlabel("time [hours]")
    ax.set_title(
        "Why the window is needed: the sample-to-sample difference collapses "
        "at every turning point",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="lower left")
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    print("wrote", out_png)


if __name__ == "__main__":
    main()
