#!/usr/bin/env python3

# Script for downloading data from an OpenScPCA data release

import argparse
import datetime
import fnmatch
import os
import pathlib
import re
import subprocess
import sys
from typing import List, Set

# enable text formatting on Windows
os.system("")

RELEASE_BUCKET = "openscpca-data-release"
RESULTS_BUCKET = "openscpca-nf-workflow-results"
TEST_BUCKET = "openscpca-temp-simdata"  # TODO: change to correct bucket


def get_releases(bucket: str, profile: str) -> List[str]:
    """
    Get the list of available releases from the OpenScPCA data bucket.
    """
    ls_cmd = ["aws", "s3", "ls", f"s3://{bucket}/"]
    if profile:
        ls_cmd += ["--profile", profile]
    if bucket == TEST_BUCKET:
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
    # get only date-based versions and remove the trailing slash
    date_re = re.compile(r"PRE\s+(20\d{2}-[01]\d-[0123]\d)/$")
    return [
        m.group(1)
        for m in (date_re.search(line) for line in ls_result.stdout.splitlines())
        if m
    ]


def get_results_modules(bucket: str, release: str, profile: str) -> List[str]:
    """
    Get the list of available results modules from the OpenScPCA results bucket.
    """
    ls_cmd = ["aws", "s3", "ls", f"s3://{bucket}/{release}/"]
    if profile:
        ls_cmd += ["--profile", profile]
    if bucket == TEST_BUCKET:
        ls_cmd += ["--no-sign-request"]
    ls_result = subprocess.run(ls_cmd, capture_output=True, text=True)
    if ls_result.returncode:
        print(
            "Error listing results modules from the OpenScPCA results bucket.",
            " Ensure you have the correct AWS permissions to access OpenScPCA data.",
            file=sys.stderr,
        )
        print(ls_result.stderr, file=sys.stderr)
        sys.exit(1)
    # get only prefixes and remove the trailing slash
    module_re = re.compile(r"PRE\s+(\w+)/")
    return [
        m.group(1)
        for m in (module_re.search(line) for line in ls_result.stdout.splitlines())
        if m
    ]


def build_sync_cmd(
    bucket: str,
    release: str,
    download_dir: pathlib.Path,
    include_patterns: List[str] = [],
    dryrun: bool = False,
    profile: str = "",
) -> List[str]:
    sync_cmd = [
        "aws",
        "s3",
        "sync",
        f"s3://{bucket}/{release}/",
        download_dir,
        "--exact-timestamps",  # replace if a file has changed at all
        "--no-progress",  # don't show progress animations
        "--exclude",
        "*",  # exclude everything by default
    ]

    for pattern in include_patterns:
        sync_cmd += ["--include", pattern]

    if dryrun:
        sync_cmd += ["--dryrun"]

    if profile:
        sync_cmd += ["--profile", profile]

    if bucket == TEST_BUCKET:
        sync_cmd += ["--no-sign-request"]

    return sync_cmd


def get_download_size(
    bucket: str,
    release: str,
    include_patterns: List[str] = ["*"],
    profile: str = "",
) -> int:
    ls_cmd = ["aws", "s3", "ls", f"s3://{bucket}/{release}/", "--recursive"]
    if profile:
        ls_cmd += ["--profile", profile]
    if bucket == TEST_BUCKET:
        ls_cmd += ["--no-sign-request"]
    file_list = subprocess.run(ls_cmd, capture_output=True, text=True)
    if file_list.returncode:
        print(
            "Error listing files in the OpenScPCA release bucket.",
            " Ensure you have the correct AWS permissions to access OpenScPCA data.",
            file=sys.stderr,
        )
        print(file_list.stderr, file=sys.stderr)
    total_size = 0
    for line in file_list.stdout.splitlines():
        print(line)
        size, file = line.split()[-2:]
        if any(fnmatch.fnmatch(file, pattern) for pattern in include_patterns):
            total_size += int(size)
    return total_size


def add_parent_dirs(patterns: List[str], dirs: List[str]) -> List[str]:
    """
    Add parent directories to each AWS include pattern.
    """
    parent_patterns = []
    for pattern in patterns:
        # Prepend only if the pattern starts with a wildcard
        if pattern.startswith("*"):
            parent_patterns += [
                # add parent directory exact or final sample in multiplexed
                f"*{d}/{pattern}"
                for d in dirs
            ]
            parent_patterns += [
                # add partial parent directory for multiplexed samples
                f"*{d}_*/{pattern}"
                for d in dirs
            ]
        else:
            parent_patterns += [pattern]
    return parent_patterns


def download_release_data(
    bucket: str,
    release: str,
    formats: Set[str],
    includes: Set[str],
    include_reports: bool,
    projects: Set[str],
    samples: Set[str],
    data_dir: pathlib.Path,
    dryrun: bool,
    profile: str,
    update_current: bool,
) -> None:
    """
    Download data for a specific release of OpenScPCA.
    """
    # Build the basic sync command
    download_dir = data_dir / release

    # Always include json, tsv metadata files and DATA_USAGE.md
    patterns = ["*.json", "*.tsv", "DATA_USAGE.md"]

    # separate bulk from other includes
    include_bulk = "bulk" in includes
    includes = includes - {"bulk"}

    if include_reports:
        patterns += ["*.html"]

    if "sce" in formats:
        patterns += [f"*_{level}.rds" for level in includes]

    if "anndata" in formats:
        patterns += [f"*_{level}_*.h5ad" for level in includes]
        patterns += [f"*_{level}_*.hdf5" for level in includes]

    # If projects or samples are specified, extend the file-only patterns to specify parent directories
    if projects:
        patterns = add_parent_dirs(patterns, projects)
    elif samples:
        patterns = add_parent_dirs(patterns, samples)

    # if samples are specified, bulk is excluded, as it is not associated with individual samples
    if include_bulk and not samples:
        if projects:
            patterns += [f"{project}/*_bulk_*.tsv" for project in projects]
        else:
            patterns += ["*_bulk_*.tsv"]

    ### build sync command and run ###
    sync_cmd = build_sync_cmd(
        bucket=bucket,
        release=release,
        download_dir=download_dir,
        include_patterns=patterns,
        dryrun=dryrun,
        profile=profile,
    )

    subprocess.run(sync_cmd, check=True)

    ### Print summary messages ###
    print("\n\n\033[1mDownload Summary\033[0m")  # bold
    print("Release:", release)
    print("Data Format:", ", ".join(formats))
    print("Processing levels:", ", ".join(includes))
    print("Projects:", ", ".join(projects) if projects else "All")
    if dryrun:
        print("Data download location:", download_dir)
    else:
        print("Downloaded data to:", download_dir)

    ### Update current link to point to new or test data if required ###
    # only do this if we are using test data or the specified release is "current" or "latest", not for specific dates
    if update_current and not dryrun:
        # update the current symlink
        current_symlink = data_dir / "current"
        current_symlink.unlink(missing_ok=True)
        current_symlink.symlink_to(release)
        print(f"Updated 'current' symlink to point to '{release}'.")


def download_results(
    bucket: str,
    release: str,
    modules: Set[str],
    data_dir: pathlib.Path,
    dryrun: bool,
    profile: str,
    update_current: bool,
) -> None:
    """
    Download results for a specific release of OpenScPCA.
    """

    download_dir = data_dir / release / "results"
    patterns = [f"{module}/*" for module in modules]

    sync_cmd = build_sync_cmd(
        bucket=bucket,
        release=release,
        download_dir=download_dir,
        include_patterns=patterns,
        dryrun=dryrun,
        profile=profile,
    )
    subprocess.run(sync_cmd, check=True)

    ### Print summary messages ###
    print("\n\n\033[1mDownload Summary\033[0m")  # bold
    print("Release:", release)
    print("Results Modules:", ", ".join(modules))
    if dryrun:
        print("Data download location:", download_dir)
    else:
        print("Downloaded data to:", download_dir)

    ### Update current link to point to new or test data  if required ###
    # only do this if we are using test data or the specified release is "current" or "latest", not for specific dates
    if update_current and not dryrun:
        # update the current symlink
        current_symlink = data_dir / "current"
        current_symlink.unlink(missing_ok=True)
        current_symlink.symlink_to(release)
        print(f"Updated 'current' symlink to point to '{release}'.")


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
        "--list-results",
        action="store_true",
        help="List the available results modules and exit.",
    )
    parser.add_argument(
        "--release",
        type=str,
        help="The release version to download. Defaults to 'current'.",
        default="",  # set in code
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
        "--module-results",
        type=str,
        default=None,
        help="The results modules to download."
        " A comma separated list of results modules to download.",
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
    anndata_formats = {"anndata", "h5ad", "hdf5"}
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
        results_bucket = TEST_BUCKET
    else:
        args.release = args.release or "current"
        bucket = RELEASE_BUCKET
        results_bucket = RESULTS_BUCKET

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

    ### List the available releases or modules ###
    all_releases = get_releases(bucket=bucket, profile=args.profile)

    # hide any future releases and sort in reverse order
    current_releases = [
        r for r in all_releases if r <= datetime.date.today().isoformat()
    ]
    current_releases.sort(reverse=True)

    # list the current releases and exit if that was what was requested
    if args.list_releases:
        print("Available release dates:\n", "\n".join(current_releases), sep="\n")
        return

    # get the release to use or exit if it is not available
    if args.test_data:
        release = "test"
    elif args.release.lower() in ["current", "latest"]:
        release = current_releases[0]
    elif args.release in all_releases:  # allow downloads from the future
        release = args.release
    else:
        print(
            f"Release dated '{args.release}' is not available. Available release dates are:\n",
            "\n".join(all_releases),
            sep="\n",
            file=sys.stderr,
        )
        sys.exit(1)

    # list results modules and exit if that was what was requested
    if args.list_results or args.module_results:
        all_modules = get_results_modules(
            bucket=results_bucket, release=release, profile=args.profile
        )
    if args.list_results:
        print(
            f"Available module results for release {release}:\n",
            "\n".join(all_modules),
            sep="\n",
        )
        return

    # check that the requested modules are available
    if args.module_results:
        modules = set(args.module_results.split(","))
        if not modules.issubset(all_modules):
            print(
                f"One or more requested modules are not available for release {release}.",
                "Available modules are:\n",
                "\n".join(all_modules),
                sep="\n",
                file=sys.stderr,
            )
            sys.exit(1)

    ### Download the data ###
    if not args.module_results:
        download_release_data(
            bucket=bucket,
            release=release,
            formats=formats,
            includes=includes,
            include_reports=args.include_reports,
            projects=args.projects.split(",") if args.projects else [],
            samples=args.samples.split(",") if args.samples else [],
            data_dir=args.data_dir,
            dryrun=args.dryrun,
            profile=args.profile,
            update_current=args.test_data
            or args.release.lower() in ["current", "latest"],
        )
    else:
        download_results(
            bucket=results_bucket,
            release=release,
            modules=modules,
            data_dir=args.data_dir,
            dryrun=args.dryrun,
            profile=args.profile,
            update_current=args.test_data
            or args.release.lower() in ["current", "latest"],
        )


if __name__ == "__main__":
    main()
