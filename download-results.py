#!/usr/bin/env python3

# Script for downloading results from the OpenScPCA-nf workflow results bucket

import argparse
import datetime
import os
import pathlib
import re
import subprocess
import sys
from typing import List, Set

# import functions from download-data.py
download_data = __import__("download-data")

# enable text formatting on Windows
os.system("")

RESULTS_BUCKET = "openscpca-nf-workflow-results"
TEST_BUCKET = "openscpca-temp-simdata"  # TODO: change to correct bucket


def get_results_modules(bucket: str, release: str, profile: str) -> List[str]:
    """
    Get the list of available results modules from the OpenScPCA workflow results bucket.
    """
    ls_cmd = ["aws", "s3", "ls", f"s3://{bucket}/{release}/"]
    if profile:
        ls_cmd += ["--profile", profile]
    if bucket == TEST_BUCKET:
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

    # get only prefixes and remove the trailing slash
    module_re = re.compile(r"PRE\s+([\S]+)/")
    return [
        m.group(1)
        for m in (module_re.search(line) for line in ls_result.stdout.splitlines())
        if m
    ]


def download_results(
    bucket: str,
    release: str,
    modules: Set[str],
    data_dir: pathlib.Path,
    projects: Set[str] = {},
    samples: Set[str] = {},
    dryrun: bool = False,
    profile: str = "",
    update_current: bool = True,
) -> None:
    """
    Download workflow results for a specific release of OpenScPCA.
    """

    download_dir = data_dir / release / "results"
    if projects and samples:
        raise ValueError("Projects and samples cannot be specified together.")

    patterns = [
        f"{m}/*{p}*" for m in modules for p in projects
    ]  # will be empty if no projects specified
    patterns += [
        f"{m}/*{s}*" for m in modules for s in samples
    ]  # will be empty if no samples specified
    if not patterns:
        patterns = [f"{m}/*" for m in modules]

    sync_cmd = download_data.build_sync_cmd(
        bucket=bucket,
        release=release,
        download_dir=download_dir,
        include_patterns=patterns,
        dryrun=dryrun,
        profile=profile,
    )
    subprocess.run(sync_cmd, check=True)

    download_size = download_data.get_download_size(
        bucket=bucket,
        release=release,
        include_patterns=patterns,
        profile=profile,
    )

    ### Print summary messages ###
    print("\n\n\033[1mDownload Summary\033[0m")  # bold
    print("Release:", release)
    print("Results Modules:", ", ".join(modules))
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
        download_data.make_size_human(download_size),
    )

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
        description="Download workflow results for OpenScPCA.",
    )
    parser.add_argument(
        "--list-releases",
        action="store_true",
        help="List the available release versions and exit.",
    )
    parser.add_argument(
        "--list-modules",
        action="store_true",
        help="List the available workflow module results and exit.",
    )
    parser.add_argument(
        "--modules",
        type=str,
        default="",
        help="The results modules to download."
        " A comma separated list of workflow modules to download results from.",
    )
    parser.add_argument(
        "--projects",
        type=str,
        default="",
        help=(
            "The project(s) to download."
            " A comma separated list of Project IDs to download results from."
            " Defaults to all. Can not be combined with `--samples`."
        ),
    )
    parser.add_argument(
        "--samples",
        type=str,
        default="",
        help="The sample(s) to download."
        " A comma separated list of Sample IDs to download results from."
        " Defaults to all. Can not be combined with `--projects`."
        " If specified, bulk files are always excluded.",
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

    args = parser.parse_args()

    ### Validate the arguments ###

    # Check that only a release or test-data was set, and set buckets and default release
    if args.release and args.test_data:
        print(
            "Only one of `--release` or `--test-data` can be set.",
            file=sys.stderr,
        )
        sys.exit(1)
    elif args.test_data:
        results_bucket = TEST_BUCKET
    else:
        args.release = args.release or "current"
        results_bucket = RESULTS_BUCKET

    # Check that projects and samples are not requested together
    if args.projects and args.samples:
        print(
            "Using both `--projects` and `--samples` options together is not supported.",
            file=sys.stderr,
        )

    # check project and sample names
    projects = {p.strip() for p in args.projects.split(",")} if args.projects else {}
    samples = {s.strip() for s in args.samples.split(",")} if args.samples else {}
    if projects and not all(p.startswith("SCPCP") for p in projects):
        print(
            "Some project ids do not start with 'SCPCP' as expected.",
            file=sys.stderr,
        )
    if samples and not all(s.startswith("SCPS") for s in samples):
        print(
            "Some sample ids do not start with 'SCPCS' as expected.",
            file=sys.stderr,
        )

    ### List the available releases or modules ###
    all_releases = download_data.get_releases(
        bucket=results_bucket, profile=args.profile
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

    # list results modules and exit if that was what was requested
    all_modules = get_results_modules(
        bucket=results_bucket, release=release, profile=args.profile
    )

    if args.list_modules:
        print(
            f"Available module results for release {release}:\n",
            "\n".join(all_modules),
            sep="\n",
        )
        return

    # check that the requested modules are available
    modules = {m.strip() for m in args.modules.split(",")} if args.modules else set()
    if not (modules and modules.issubset(all_modules)):
        if modules:
            print(
                f"One or more requested modules are not available for release {release}.",
                file=sys.stderr,
            )
        else:
            print(
                "At least one results module must be specified with --modules.",
                file=sys.stderr,
            )
        print(
            "\nAvailable modules are:",
            "\n".join(all_modules),
            sep="\n",
            file=sys.stderr,
        )
        sys.exit(1)

    ### Download the results files ###
    download_results(
        bucket=results_bucket,
        release=release,
        modules=modules,
        data_dir=args.data_dir,
        projects=projects,
        samples=samples,
        dryrun=args.dryrun,
        profile=args.profile,
        update_current=args.test_data or args.release.lower() in ["current", "latest"],
    )


if __name__ == "__main__":
    main()
