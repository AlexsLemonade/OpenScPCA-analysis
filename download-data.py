#!/usr/bin/env python3

# Script for downloading data from an OpenScPCA data release

import argparse
import os
import pathlib
import re
import subprocess
import sys

# enable text formatting on Windows
os.system("")


def add_parent_dirs(patterns: list[str], dirs: list[str]) -> list[str]:
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
        "--bucket",
        type=str,
        help="The S3 bucket to download the data from. Default is OpenScPCA data release bucket.",
        default="analysis-s3-992382809252-us-east-2",  # TODO: change to correct bucket
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress output except errors.",
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
    ls_cmd = ["aws", "s3", "ls", f"s3://{args.bucket}/"]
    ls_result = subprocess.run(ls_cmd, capture_output=True, text=True)
    if ls_result.returncode:
        print(
            f"Error listing releases from S3 bucket '{args.bucket}'.",
            "Ensure you have the correct permissions and the bucket exists.",
            file=sys.stderr,
        )
        print(ls_result.stderr, file=sys.stderr)
        sys.exit(1)
    # get only date-based versions and remove the trailing slash
    date_re = re.compile(r"(\d{4}-\d{2}-\d{2})/?")
    releases = [
        m.group(1)
        for m in (date_re.search(line) for line in ls_result.stdout.splitlines())
        if m
    ]

    # list the releases and exit if that was what was requested
    if args.list_releases:
        print("Available release dates:", "\n".join(releases), sep="\n")
        return

    # get the release to use or exit if it is not available
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

    ### Download the data ###
    # Build the basic sync command
    sync_cmd = [
        "aws",
        "s3",
        "sync",
        f"s3://{args.bucket}/{release}/",
        f"{args.data_dir}/{release}/",
        "--exact-timestamps",  # replace if a file has changed at all
        "--exclude",
        "*",  # exclude everything by default
    ]

    # Always include report files
    patterns = ["*.html", "*.json"]

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

    if args.quiet:
        sync_cmd += ["--only-show-errors"]

    subprocess.run(sync_cmd, check=True)

    ### Print summary messages ###
    if not args.quiet:
        print("\033[1mDownload Summary\033[0m")  # bold
        print("Release:", release)
        print("Data Format:", ", ".join(formats))
        print("Processing levels:", ", ".join(includes))
        print("Downloaded data to:", args.data_dir / release)

    ### Update current link to point to the new data release ###
    # only do this if the specified release is "current" or "latest", not for specific dates
    if args.release.lower() in ["current", "latest"] and not args.dryrun:
        # update the current symlink
        current_symlink = args.data_dir / "current"
        current_symlink.unlink(missing_ok=True)
        current_symlink.symlink_to(args.data_dir / release)
        if not args.quiet:
            print(f"Updated 'current' symlink to point to '{release}'.")


if __name__ == "__main__":
    main()
