import numpy as np
import argparse
from pathlib import Path
import pandas as pd
import sys

AMR_WIND_PATH = Path(__file__).parents[2].resolve()
sys.path.append(str(AMR_WIND_PATH / "tools"))
import amrex_particle
from amrex_particle import AmrexParticleFile


def main():
    """Compare particle files"""
    parser = argparse.ArgumentParser(description="A comparison tool for particles")
    parser.add_argument("f0", help="A particle directory", type=str)
    parser.add_argument("f1", help="A particle directory", type=str)
    parser.add_argument(
        "-a",
        "--abs_tol",
        help="Absolute tolerance (default is 0)",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "-r",
        "--rel_tol",
        help="Relative tolerance (default is 0)",
        type=float,
        default=0.0,
    )
    args = parser.parse_args()

    assert Path(args.f0).is_dir()
    assert Path(args.f1).is_dir()

    cols_to_drop = ["uid", "set_id", "probe_id"]
    p0df = AmrexParticleFile(args.f0)().drop(cols_to_drop, axis=1)
    p1df = AmrexParticleFile(args.f1)().drop(cols_to_drop, axis=1)

    adiff = np.sqrt(np.square(p0df - p1df).sum(axis=0))
    rdiff = np.sqrt(np.square(p0df - p1df).sum(axis=0)) / np.sqrt(
        np.square(p0df).sum(axis=0)
    )
    adiff = adiff.to_frame(name="absolute_error")
    rdiff = rdiff.to_frame(name="relative_error")
    diff = pd.concat([adiff, rdiff], axis=1).fillna(value={"relative_error": 0.0})

    print(diff)
    assert (diff["absolute_error"] <= args.abs_tol).all()
    assert (diff["relative_error"] <= args.rel_tol).all()


if __name__ == "__main__":
    main()
