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
TEST_BUCKET = "openscpca-test-data-release-public-access"


def get_releases(bucket: str, profile: str, test_data: bool) -> List[str]:
    """
    Get the list of available releases from an OpenScPCA bucket.
    """
    ls_cmd = ["aws", "s3", "ls", f"s3://{bucket}/"]
    if profile:
        ls_cmd += ["--profile", profile]
    if test_data:
        ls_cmd += ["--no-sign-request"]
    ls_result = subprocess.run(ls_cmd, capture_output=True, text=True)
    if ls_result.returncode:  # authentication errors, usually
        print(
            "Error listing release versions from the OpenScPCA bucket.\n\n"
            "Make sure you have the correct profile active (or use the --profile option)"
            " and run `aws sso login` before running this script.\n\n"
            "AWS Error: ",
            ls_result.stderr,  # print the AWS error text too
            file=sys.stderr,
        )
        sys.exit(ls_result.returncode)

    # get only date-based versions and remove the trailing slash
    date_re = re.compile(r"PRE\s+(20\d{2}-[01]\d-[0123]\d)/$")
    return [
        m.group(1)
        for m in (date_re.search(line) for line in ls_result.stdout.splitlines())
        if m
    ]


def build_sync_cmd(
    bucket: str,
    release: str,
    test_data: bool,
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
    ]

    if include_patterns:
        sync_cmd += ["--exclude", "*"]

    for pattern in include_patterns:
        sync_cmd += ["--include", pattern]

    if dryrun:
        sync_cmd += ["--dryrun"]

    if profile:
        sync_cmd += ["--profile", profile]

    if test_data:
        sync_cmd += ["--no-sign-request"]

    return sync_cmd


def get_download_size(
    bucket: str,
    release: str,
    test_data: bool,
    include_patterns: List[str] = ["*"],
    profile: str = "",
) -> int:
    """
    Get the total size of files that will be downloaded from AWS S3 that match the include patterns.
    """
    ls_cmd = ["aws", "s3", "ls", f"s3://{bucket}/{release}/", "--recursive"]
    if profile:
        ls_cmd += ["--profile", profile]
    if test_data:
        ls_cmd += ["--no-sign-request"]
    file_list = subprocess.run(ls_cmd, capture_output=True, text=True)
    file_list.check_returncode()

    total_size = 0
    for line in file_list.stdout.splitlines():
        size, file = line.split()[-2:]
        if any(
            fnmatch.fnmatch(file, release + "/" + pattern)
            for pattern in include_patterns
        ):
            total_size += int(size)
    return total_size


def make_size_human(size: int) -> str:
    """
    Convert a size in bytes to something human readable.
    """
    for unit in ["B", "KiB", "MiB", "GiB", "TiB"]:
        if size < 1024.0 or unit == "TiB":
            break
        size /= 1024.0
    if unit == "B":
        return f"{size:.0f} B"
    return f"{size:.2f} {unit}"


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


def update_symlink(
    data_dir: pathlib.Path,
    target: str
) -> None:
    """
    Update symlink to direct to target
    """
    current_symlink = data_dir / "current"
    current_symlink.unlink(missing_ok=True)
    current_symlink.symlink_to(target)
    print(f"Updated 'current' symlink to point to '{target}'.")



def download_release_data(
    bucket: str,
    release: str,
    test_data: bool,
    data_dir: pathlib.Path,
    formats: Set[str],
    stages: Set[str],
    include_reports: bool = False,
    projects: Set[str] = {},
    samples: Set[str] = {},
    dryrun: bool = False,
    profile: str = "",
    update_current: bool = True,
) -> None:
    """
    Download data for a specific release of OpenScPCA.
    """
    # Build the basic sync command
    download_dir = data_dir / release

    # Always include json, tsv metadata files and DATA_USAGE.md
    patterns = ["*.json", "*_metadata.tsv", "DATA_USAGE.md"]

    # separate bulk from other stages
    include_bulk = "bulk" in stages
    stages = stages - {"bulk"}

    if include_reports:
        patterns += ["*.html"]

    if "sce" in formats:
        patterns += [f"*_{level}.rds" for level in stages]

    if "anndata" in formats:
        patterns += [f"*_{level}_*.h5ad" for level in stages]
        patterns += [f"*_{level}_*.hdf5" for level in stages]

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
        test_data=test_data,
        download_dir=download_dir,
        include_patterns=patterns,
        dryrun=dryrun,
        profile=profile,
    )

    subprocess.run(sync_cmd, check=True)

    download_size = get_download_size(
        bucket=bucket,
        release=release,
        test_data=test_data,
        include_patterns=patterns,
        profile=profile,
    )

    ### Print summary messages ###
    print("\n\n\033[1mDownload Summary\033[0m")  # bold
    print("Release:", release)
    if formats:
        print("Data Format:", ", ".join(formats))
    if stages:
        print("Processing levels:", ", ".join(stages))
    if projects:
        print("Projects:", ", ".join(projects))
    if samples:
        print("Samples:", ", ".join(samples))
    if dryrun:
        print("Data download location:", download_dir)
    else:
        print("Downloaded data to:", download_dir)
    print(
        "Total download size (includes previous downloads):",
        make_size_human(download_size),
    )

    ### Update current link to point to new or test data if required ###
    # only do this if we are using test data or the specified release is "current" or "latest", not for specific dates
    if update_current:
        if not dryrun:
            update_symlink(download_dir, release)
        else:
            print(f"\nIf this were not a dry run, the 'current' symlink would be updated to point to '{release}'.")



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
        "--projects",
        type=str,
        default="",
        help=(
            "The project(s) to download."
            " A comma separated list of Project IDs to download."
            " Defaults to all. Can not be combined with `--samples`."
        ),
    )
    parser.add_argument(
        "--samples",
        type=str,
        default="",
        help="The sample(s) to download."
        " A comma separated list of Sample IDs to download."
        " Defaults to all. Can not be combined with `--projects`."
        " If specified, bulk files are always excluded.",
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
        "--process-stage",
        type=str,
        default="processed",
        help=(
            "The stage of processing for the data files to include."
            " One or more of 'processed', 'filtered', 'unfiltered', or 'bulk'."
            " Defaults to 'processed'."
            " For more than one level, use a comma separated list."
        ),
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
        " To switch back, rerun this script with either the `--release current` option, or the `--update-symlink release` option.",
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Download only the metadata files and not the data files."
        " To also download QC reports, combine with the --include-reports option."
        " Can be combined with --projects, but not with --samples. --format and --process-stage are ignored.",
    )
    parser.add_argument(
        "--include-reports",
        action="store_true",
        help="Include html report files in the download."
        " Note that test data does not include report files.",
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
        "--dryrun",
        action="store_true",
        help="Perform a dry run of the download: show what would be done but do not download anything.",
    )
    parser.add_argument(
        "--profile",
        type=str,
        default="",
        help="The AWS profile to use for the download. Uses the current default profile if undefined.",
    )
    parser.add_argument(
        "--update-symlink",
        type=str,
        default="",
        help="The release to update the 'current' symlink to direct to. Data will not be re-downloaded data if it is already present."
    )
    args = parser.parse_args()

    ### Validate the arguments ###
    validation_error = False
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
        validation_error = True

    # Check include levels are valid & make a set
    process_stages = {"unfiltered", "filtered", "processed", "bulk"}
    stages = {x.strip().lower() for x in args.process_stage.split(",")}
    if not all(x in process_stages for x in stages):
        print(
            f"process-stage option '{args.process_stage}' is not valid.",
            "Must be 'processed', 'filtered', 'unfiltered', 'bulk', or a comma separated list of those.",
            file=sys.stderr,
        )
        validation_error = True

    # Check that only a release or test-data was set, and set buckets and default release
    if args.release and args.test_data:
        print(
            "Only one of `--release` or `--test-data` can be set.",
            file=sys.stderr,
        )
        validation_error = True
    elif args.test_data:
        bucket = TEST_BUCKET
    else:
        args.release = args.release or "current"
        bucket = RELEASE_BUCKET
    # Check that projects and samples are not requested together
    if args.projects and args.samples:
        print(
            "Using both `--projects` and `--samples` options together is not supported.",
            file=sys.stderr,
        )
        validation_error = True

    if args.metadata_only and args.samples:
        print(
            "Using `--metadata-only` and `--samples` options together is not supported.",
            file=sys.stderr,
        )
        validation_error = True

    if validation_error:
        sys.exit(1)

    # check project and sample names
    projects = {p.strip() for p in args.projects.split(",")} if args.projects else {}
    samples = {s.strip() for s in args.samples.split(",")} if args.samples else {}
    if projects and not all(p.startswith("SCPCP") for p in projects):
        print(
            "Some project ids do not start with 'SCPCP' as expected.",
            file=sys.stderr,
        )
    if samples and not all(s.startswith("SCPCS") for s in samples):
        print(
            "Some sample ids do not start with 'SCPCS' as expected.",
            file=sys.stderr,
        )

    if args.samples and "bulk" in stages:
        print(
            "Bulk data is not available for individual samples, so bulk data will be skipped.",
            file=sys.stderr,
        )


    ### Update symlink if requested, without downloading data, if the target directory exists
    if args.update_symlink:
        current_symlink = args.data_dir / "current"
        desired_target = args.data_dir / args.update_symlink
        if desired_target.exists():
            update_symlink(args.data_dir, args.update_symlink)
            sys.exit(0)
        else:
            print(
                f"The requested target directory {args.update_symlink} for the 'current' symlink does not exist.\n",
                " Please instead download the desired release with the `--release` flag.",
                file=sys.stderr,
            )
            sys.exit(1)

    ### List the available releases or modules ###
    all_releases = get_releases(
        bucket=bucket, profile=args.profile, test_data=args.test_data
    )

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

    if args.metadata_only:
        formats = set()
        stages = set()

    ### Download the data ###
    download_release_data(
        bucket=bucket,
        release=release,
        test_data=args.test_data,
        data_dir=args.data_dir,
        formats=formats,
        stages=stages,
        include_reports=args.include_reports,
        projects=args.projects.split(",") if args.projects else [],
        samples=args.samples.split(",") if args.samples else [],
        dryrun=args.dryrun,
        profile=args.profile,
        update_current=args.test_data or args.release.lower() in ["current", "latest"],
    )


if __name__ == "__main__":
    main()
