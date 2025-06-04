# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "matplotlib",
#     "numpy",
#     "scipy",
#     "rasterio",
# ]
# ///

# Usage: uv run --script unit_tests/utilities/test_raster_asc.py --plot
import argparse

import matplotlib.pyplot as plt
import numpy as np
import rasterio
from scipy.interpolate import RegularGridInterpolator


def check_rasterio_cell_values(raster_path, points):
    with rasterio.open(raster_path) as src:
        coords = [(x, y) for x, y, _ in points]
        for idx, (val, (x, y, expected)) in enumerate(zip(src.sample(coords), points)):
            v = val[0]
            print(
                f"Point [{idx}] (x={x:.3f}, y={y:.3f}) -> rasterio.sample = {v:.8f}, expected = {expected:.8f}"
            )


def main(plot=False, npoints=20, seed=42):
    # --- 1. Generate raster with Gaussian hill + asymmetry ---
    ncols = nrows = 50
    xllcorner = yllcorner = 0.0
    cellsize = 1.0
    nodata = -9999.0

    x = np.arange(ncols)
    y = np.arange(nrows)
    X, Y = np.meshgrid(x, y)

    # Gaussian hill parameters
    x0 = y0 = (ncols - 1) / 2.0
    sigma = 10.0
    hill = np.exp(-((X - x0) ** 2 + (Y - y0) ** 2) / (2 * sigma**2))

    # Add asymmetry: linear gradient in x and y
    asym = 0.5 * (X / ncols) + 0.2 * (Y / nrows)
    values = hill + asym

    # --- 2. Write to raster.asc (ArcGIS ASCII format) ---
    with open("unit_tests/utilities/raster.asc", "w") as f:
        f.write(f"ncols         {ncols}\n")
        f.write(f"nrows         {nrows}\n")
        f.write(f"xllcorner     {xllcorner}\n")
        f.write(f"yllcorner     {yllcorner}\n")
        f.write(f"cellsize      {cellsize}\n")
        f.write(f"NODATA_value  {nodata}\n")
        for row in values[::-1]:
            f.write(" ".join(f"{val:.6f}" for val in row) + "\n")

    # --- 3. Read the raster.asc back ---
    def load_ascii_raster(filename):
        with open(filename) as f:
            header = {}
            for _ in range(6):
                key, val = f.readline().split()
                header[key.lower()] = (
                    float(val) if "." in val or "e" in val.lower() else int(val)
                )
            data = np.loadtxt(f)
        arr = np.flipud(data)
        return arr, header

    arr, header = load_ascii_raster("unit_tests/utilities/raster.asc")
    ncols, nrows = header["ncols"], header["nrows"]
    xll, yll, dx = header["xllcorner"], header["yllcorner"], header["cellsize"]

    # --- 4. Interpolate at random points and plot ---
    x = np.linspace(xll, xll + dx * (ncols - 1), ncols)
    y = np.linspace(yll, yll + dx * (nrows - 1), nrows)
    interp = RegularGridInterpolator((y, x), arr)

    np.random.seed(seed)
    rand_pts = np.random.rand(npoints, 2)
    rand_pts[:, 0] = rand_pts[:, 0] * (x[-1] - x[0]) + x[0]
    rand_pts[:, 1] = rand_pts[:, 1] * (y[-1] - y[0]) + y[0]
    rand_vals = interp(rand_pts[:, [1, 0]])  # (y, x) order

    if plot:
        plt.figure(figsize=(7, 6))
        plt.imshow(
            arr,
            origin="lower",
            extent=(xll, xll + dx * ncols, yll, yll + dx * nrows),
            cmap="viridis",
        )
        plt.colorbar(label="Value")
        plt.title("Asymmetric Gaussian Raster")
        plt.scatter(
            rand_pts[:, 0],
            rand_pts[:, 1],
            c=rand_vals,
            edgecolor="k",
            cmap="cool",
            s=80,
            label="Random points",
        )
        # Annotate each point with its index
        for i, (xpt, ypt) in enumerate(zip(rand_pts[:, 0], rand_pts[:, 1])):
            plt.annotate(
                str(i),
                (xpt, ypt),
                textcoords="offset points",
                xytext=(5, 5),
                fontsize=9,
                color="white",
                weight="bold",
            )
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.savefig("raster_plot.png")
        plt.show()

    # --- Print for C++ ---
    sample_points = list(zip(rand_pts[:, 0], rand_pts[:, 1], rand_vals))
    print("std::vector<std::tuple<double, double, double>> sample_points = {")
    for xi, yi, vi in sample_points:
        print(f"    {{{xi:.8f}, {yi:.8f}, {vi:.8f}}},")
    print("};")

    check_rasterio_cell_values("unit_tests/utilities/raster.asc", sample_points)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot", action="store_true", help="Show plot")
    parser.add_argument(
        "--npoints", type=int, default=20, help="Number of random points"
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()
    main(plot=args.plot, npoints=args.npoints, seed=args.seed)
