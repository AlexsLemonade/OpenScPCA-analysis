#!/usr/bin/env python3

# Script for downloading data from an OpenScPCA data release

import argparse
import pathlib
import re
import subprocess

# TODO: change to correct bucket
aws_bucket = "analysis-s3-992382809252-us-east-2"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download data for OpenScPCA.",
    )
    parser.add_argument(
        "--list-releases",
        action="store_true",
        help="List the available release versions.",
    )
    parser.add_argument(
        "--release",
        type=str,
        help="The release version to download.",
        default="current",
    )
    parser.add_argument(
        "--format",
        type=str,
        help=(
            "The format to download the data in. Either 'SCE' or 'AnnData'. Defaults to 'SCE'."
            "For more than one format, use a comma separated list with no spaces."
        ),
        default="SCE",
    )
    parser.add_argument(
        "--data-dir",
        type=pathlib.Path,
        help="The directory to download the data to. Defaults to the `data` directory in this repository.",
        default=pathlib.Path(__file__).parent / "data",
    )

    args = parser.parse_args()

    # Check formats are valid
    sce_formats = ["sce", "rds"]
    anndata_formats = ["anndata", "hdf5", "h5ad"]
    formats = [f.lower() for f in args.format.split(",")]
    if not all(f in sce_formats + anndata_formats for f in formats):
        raise ValueError(
            f"Format '{args.format}' not available. Must be 'SCE', 'AnnData', or a comma separated list."
        )

    aws_get_prefixes = f"aws s3 ls 's3://{aws_bucket}/'" "| awk '{print $NF}'"
    prefixes = subprocess.run(
        aws_get_prefixes, shell=True, capture_output=True, text=True
    ).stdout.split()
    # get only date-based versions and remove the trailing slash
    releases = [p.strip("/") for p in prefixes if re.match(r"\d{4}-\d{2}-\d{2}", p)]

    if args.list_releases:
        print("Available release dates:")
        print("\n".join(releases))
        return

    # get the release to use
    if args.release.lower() in ["current", "latest"]:
        release = max(releases)
    elif args.release in releases:
        release = args.release
    else:
        raise ValueError(
            f"Release dated '{args.release}' not available. Available dates are: {releases}."
        )

    # download the data
    # always include the metadata
    aws_sync_command = [
        "aws",
        "s3",
        "sync",
        f"s3://{aws_bucket}/{release}/",
        f"{args.data_dir}/{release}/",
        "--exact-timestamps",  # replace if the files has changed
        "--exclude",
        "*",  # exclude everything
        "--include",
        "*.html",  # always include html reports
        "--include",
        "*.json",  # always include json metadata
    ]

    subprocess.run(aws_sync_command)


if __name__ == "__main__":
    main()
