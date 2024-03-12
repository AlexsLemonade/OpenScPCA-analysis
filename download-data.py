#!/usr/bin/env python3

# Script for downloading data from an OpenScPCA data release

import argparse
import pathlib
import re
import subprocess
import sys

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
        help="The release version to download. Defaults to 'current'.",
        default="current",
    )
    parser.add_argument(
        "--format",
        type=str,
        help=(
            "The format to download the data in. Either 'SCE' or 'AnnData'."
            " Defaults to 'SCE'."
            " For more than one format, use a comma separated list with no spaces."
        ),
        default="SCE",
    )
    parser.add_argument(
        "--include",
        type=str,
        default="processed",
        help=(
            "The level of processing for the experiment files to include."
            " One or more of 'processed', 'filtered', 'unfiltered', or 'bulk'."
            " Defaults to 'processed'."
            " For more than one level, use a comma separated list with no spaces."
        ),
    )
    parser.add_argument(
        "--projects",
        type=str,
        default=None,
        help=(
            "The project(s) to download."
            " A comma separated list of Project IDs to download."
            " Defaults to all. Can not be combined with `--samples`."
        ),
    )
    parser.add_argument(
        "--samples",
        type=str,
        default=None,
        help="The sample(s) to download."
        " A comma separated list of Sample IDs to download."
        " Defaults to all. Can not be combined with `--projects`."
        " If specified, bulk files are always excluded.",
    )
    parser.add_argument(
        "--data-dir",
        type=pathlib.Path,
        help=(
            "The directory to download the data to."
            " Defaults to the `data` directory in this repository."
        ),
        default=pathlib.Path(__file__).parent / "data",
    )

    args = parser.parse_args()

    # Check formats are valid
    sce_formats = {"sce", "rds"}
    anndata_formats = {"anndata", "hdf5", "h5ad"}
    formats = {f.lower() for f in args.format.split(",")}
    if not all(f in sce_formats | anndata_formats for f in formats):
        print(
            f"Format '{args.format}' not available.",
            "Must be 'SCE', 'AnnData', or a comma separated list of those options.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Check include levels are valid
    include_levels = {"unfiltered", "filtered", "processed", "bulk"}
    includes = {x.lower() for x in args.include.split(",")}
    if not all(x in include_levels for x in includes):
        print(
            f"Include option '{args.include}' is not valid.",
            "Must be 'processed', 'filtered','unfiltered', 'filtered', 'bulk', or a comma separated list of those.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Check that projects and samples are not requested together
    if args.projects and args.samples:
        print(
            "Using both `--projects` and `--samples` options together is not supported.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.samples and "bulk" in includes:
        print(
            "Bulk data is not available for individual samples, so bulk data will be skipped.",
            file=sys.stderr,
        )

    # pull bulk to its own variable and remove from includes
    include_bulk = "bulk" in includes
    includes = includes - {"bulk"}

    ls_cmd = ["aws", "s3", "ls", f"s3://{aws_bucket}/"]
    ls_result = subprocess.run(ls_cmd, capture_output=True, text=True)
    ls_result.check_returncode()
    # get only date-based versions and remove the trailing slash
    date_re = re.compile(r"(\d{4}-\d{2}-\d{2})/?")
    releases = [
        m.group(1)
        for m in (date_re.search(line) for line in ls_result.stdout.splitlines())
        if m
    ]

    if args.list_releases:
        print(
            "Available release dates:", "\n".join(releases), sep="\n", file=sys.stderr
        )
        return

    # get the release to use
    if args.release.lower() in ["current", "latest"]:
        release = max(releases)
    elif args.release in releases:
        release = args.release
    else:
        print(
            f"Release dated '{args.release}' is not available. Available release dates are:",
            "\n".join(releases),
            sep="\n",
            file=sys.stderr,
        )
        sys.exit(1)

    # download the data
    # always include the metadata
    sync_cmd = [
        "aws",
        "s3",
        "sync",
        f"s3://{aws_bucket}/{release}/",
        f"{args.data_dir}/{release}/",
        "--exact-timestamps",  # replace if the files has changed
        "--exclude",
        "*",  # exclude everything by default
    ]

    # Always include report files
    patterns = ["*.html", "*.json"]

    if "sce" in formats:
        patterns.extend([f"*_{level}.rds" for level in includes])

    if "anndata" in formats:
        patterns.extend([f"*_{level}_*.hdf5" for level in includes])

    if args.projects:
        project_patterns = []
        for project in args.projects.split(","):
            project_patterns.extend([f"*/{project}/*/{p}" for p in patterns])
        patterns = project_patterns

    if args.samples:
        sample_patterns = []
        for sample in args.samples.split(","):
            sample_patterns.extend([f"*/{sample}/{p}" for p in patterns])
        patterns = sample_patterns

    if include_bulk and not args.samples:
        # if samples are specified, bulk is excluded, as it is not associated with samples
        if args.projects:
            patterns.extend(
                [f"*/{project}/*_bulk_*.tsv" for project in args.projects.split(",")]
            )
        else:
            patterns.extend(["*_bulk_*.tsv"])

    for p in patterns:
        sync_cmd.extend(["--include", p])

    sync_cmd.append("--dryrun")
    print(" ".join(sync_cmd))
    ## TODO: add report about what was done


if __name__ == "__main__":
    main()
