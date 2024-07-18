#!/usr/bin/env python3

"""
This script runs the SEACell algorithm on the given dataset.
"""

import argparse
import pathlib

import AnnData
import SEACells


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the SEACell algorithm on the given dataset."
    )
    parser.add_argument(
        "datafile", type=pathlib.Path, help="The input data in H5AD format."
    )
    parser.add_argument("output", type=pathlib.Path, help="The output filename.")

    args = parser.parse_args()

    adata = AnnData.read(args.datafile)


if __name__ == "__main__":
    main()
