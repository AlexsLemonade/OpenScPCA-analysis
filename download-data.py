#!/usr/bin/env python3

# Script for downloading data from an OpenScPCA data release

import argparse
import datetime
import os
import pathlib
import re
import subprocess
import sys
from typing import List

# enable text formatting on Windows
os.system("")

RELEASE_BUCKET = "openscpca-data-release"
TEST_BUCKET = "openscpca-temp-simdata"  # TODO: change to correct bucket


def add_parent_dirs(patterns: List[str], dirs: List[str]) -> List[str]:
    """
    Add parent directories to each AWS include pattern.
    """
    parent_patterns = []
    for pattern in patterns:
        # Prepend only if the pattern starts with a wildcard
        if pattern.startswith("*"):
            parent_patterns += [f"*{d}/{pattern}" for d in dirs]
        else:
            parent_patterns += [pattern]
    return parent_patterns


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download data for OpenScPCA.",
    )
    parser.add_argument(
        "--list-releases",
        action="store_true",
        help="List the available release versions and exit.",
    )
    parser.add_argument(
        "--release",
        type=str,
        help="The release version to download. Defaults to 'current'.",
        default=None,  # default is set in the script
    )
    parser.add_argument(
        "--test-data",
        action="store_true",
        help="Download test data from the test bucket and direct the `current` symlink to the test data directory."
        " To switch back, rerun this script with the `--release current` option.",
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
        "--include-reports",
        action="store_true",
        help="Include html report files in the download.",
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
        "--dryrun",
        action="store_true",
        help="Perform a dry run of the download: show what would be done but do not download anything.",
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
    parser.add_argument(
        "--profile",
        type=str,
        default=None,
        help="The AWS profile to use for the download. Uses the current default profile if undefined.",
    )

    args = parser.parse_args()

    ### Validate the arguments ###

    # Check formats are valid and make a set
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

    # Check include levels are valid & make a set
    include_levels = {"unfiltered", "filtered", "processed", "bulk"}
    includes = {x.lower() for x in args.include.split(",")}
    if not all(x in include_levels for x in includes):
        print(
            f"Include option '{args.include}' is not valid.",
            "Must be 'processed', 'filtered','unfiltered', 'filtered', 'bulk', or a comma separated list of those.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Check that only a release to test-data was set, and set buckets and default release
    if args.release and args.test_data:
        print(
            "Only one of `--release` or `--test-data` can be set.",
            file=sys.stderr,
        )
        sys.exit(1)
    elif args.test_data:
        bucket = TEST_BUCKET
    else:
        bucket = RELEASE_BUCKET
        args.release = args.release or "current"

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

    # pull bulk to its own variable and remove from includes set
    include_bulk = "bulk" in includes
    includes = includes - {"bulk"}

    ### List the available releases ###
    ls_cmd = ["aws", "s3", "ls", f"s3://{bucket}/"]
    if args.profile:
        ls_cmd += ["--profile", args.profile]
    if args.test_data:
        ls_cmd += ["--no-sign-request"]
    ls_result = subprocess.run(ls_cmd, capture_output=True, text=True)
    if ls_result.returncode:
        print(
            "Error listing release versions from the OpenScPCA release bucket.",
            " Ensure you have the correct AWS permissions to access OpenScPCA data.",
            file=sys.stderr,
        )
        print(ls_result.stderr, file=sys.stderr)
        sys.exit(1)
    # get only date-based versions or "test" and remove the trailing slash
    date_re = re.compile(r"((\d{4}-\d{2}-\d{2})|(test))/?")
    all_releases = [
        m.group(1)
        for m in (date_re.search(line) for line in ls_result.stdout.splitlines())
        if m
    ]

    # hide any future releases
    current_releases = [r for r in all_releases if r <= datetime.date.today().isoformat()]

    # list the current releases and exit if that was what was requested
    if args.list_releases:
        print("Available release dates:", "\n".join(current_releases), sep="\n")
        return

    # get the release to use or exit if it is not available
    if args.test_data:
        release = "test"
    elif args.release.lower() in ["current", "latest"]:
        release = max(current_releases)
    elif args.release in all_releases:  # allow downloads from the future
        release = args.release
    else:
        print(
            f"Release dated '{args.release}' is not available. Available release dates are:",
            "\n".join(all_releases),
            sep="\n",
            file=sys.stderr,
        )
        sys.exit(1)

    ### Download the data ###
    # Build the basic sync command
    sync_cmd = [
        "aws",
        "s3",
        "sync",
        f"s3://{bucket}/{release}/",
        f"{args.data_dir}/{release}/",
        "--exact-timestamps",  # replace if a file has changed at all
        "--no-progress",  # don't show progress animations
        "--exclude",
        "*",  # exclude everything by default
    ]

    # Always include json, tsv metadata files and DATA_USAGE.md
    patterns = ["*.json", "*.tsv", "DATA_USAGE.md"]

    if args.include_reports:
        patterns += ["*.html"]

    if "sce" in formats:
        patterns += [f"*_{level}.rds" for level in includes]

    if "anndata" in formats:
        patterns += [f"*_{level}_*.hdf5" for level in includes]

    # If projects or samples are specified, extend the file-only patterns to specify parent directories
    if args.projects:
        patterns = add_parent_dirs(patterns, args.projects.split(","))
    elif args.samples:
        patterns = add_parent_dirs(patterns, args.samples.split(","))

    # if samples are specified, bulk is excluded, as it is not associated with individual samples
    if include_bulk and not args.samples:
        if args.projects:
            patterns += [
                f"{project}/*_bulk_*.tsv" for project in args.projects.split(",")
            ]
        else:
            patterns += ["*_bulk_*.tsv"]

    ### Add patterns to the sync command and run it! ###
    for p in patterns:
        sync_cmd += ["--include", p]

    if args.dryrun:
        sync_cmd += ["--dryrun"]

    if args.profile:
        sync_cmd += ["--profile", args.profile]

    if args.test_data:
        sync_cmd += ["--no-sign-request"]

    subprocess.run(sync_cmd, check=True)

    ### Print summary messages ###
    print("\n\n\033[1mDownload Summary\033[0m")  # bold
    print("Release:", release)
    print("Data Format:", ", ".join(formats))
    print("Processing levels:", ", ".join(includes))
    print("Downloaded data to:", args.data_dir / release)

    ### Update current link to point to new or test data ###
    # only do this if we are using test data or the specified release is "current" or "latest", not for specific dates
    if (
        args.test_data or args.release.lower() in ["current", "latest"]
    ) and not args.dryrun:
        # update the current symlink
        current_symlink = args.data_dir / "current"
        current_symlink.unlink(missing_ok=True)
        current_symlink.symlink_to(release)
        print(f"Updated 'current' symlink to point to '{release}'.")


if __name__ == "__main__":
    main()
