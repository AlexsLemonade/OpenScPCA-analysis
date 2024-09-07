#!/usr/bin/env python3

# Script for syncing data from an OpenScPCA analysis module to an S3 bucket

import argparse
import os
import pathlib
import subprocess
import sys

# enable text formatting on Windows
os.system("")


def find_git_root() -> pathlib.Path:
    """Find the root directory of the git repository from the calling location"""
    repo_root = pathlib.Path.cwd()
    while not (repo_root / ".git").is_dir():  # search for the .git directory
        repo_root = repo_root.parent
        if repo_root == "/":
            raise FileNotFoundError("Could not find the repository root directory")
    return repo_root


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Sync data from an analysis module to your S3 bucket.",
    )
    parser.add_argument(
        "module",
        metavar="MODULE",
        type=str,
        help="The name of the analysis module to sync results from.",
    )
    parser.add_argument(
        "--bucket",
        "-b",
        type=str,
        required=False,
        help="The name of the S3 bucket to sync results to. Will use the OPENSCPCA_RESULTS_BUCKET environment variable if not specified.",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_false",
        dest="sync_plots",
        help="Skip syncing the plots directory.",
    )
    parser.add_argument(
        "--skip-results",
        action="store_false",
        dest="sync_results",
        help="Skip syncing the results directory.",
    )
    parser.add_argument(
        "--destructive-sync",
        action="store_true",
        dest="delete_extra",
        help=(
            "Delete any files in the S3 bucket that are not present in the local directories being synced."
            " Default behavior is to leave any additional files in the S3 bucket."
        ),
    )
    parser.add_argument(
        "--dryrun",
        action="store_true",
        help="Perform a dry run: show what would be done but do not copy or delete any files.",
    )
    parser.add_argument(
        "--profile",
        type=str,
        default=None,
        help="The AWS profile to use. Uses the current default profile if undefined.",
    )

    args = parser.parse_args()

    # set the bucket name (if not specified, use an environment variable)
    if args.bucket:
        bucket = args.bucket
    else:  # if no bucket is specified, use the environment
        bucket = os.environ.get("OPENSCPCA_RESULTS_BUCKET")
    if not bucket:
        print(
            "No bucket specified and OPENSCPCA_RESULTS_BUCKET environment variable is not set.",
            file=sys.stderr,
        )
        sys.exit(1)

    ## Set up the paths and directories ##
    # Find the repository root directory
    repo_root = find_git_root()

    # set module path (using pathlib)
    module_root = repo_root / "analyses" / args.module

    # check that module exists
    if not module_root.exists():
        print(
            f"Analysis module '{args.module}' does not exist within the 'analyses' directory.",
            file=sys.stderr,
        )
        sys.exit(1)

    sync_cmd = [
        "aws",
        "s3",
        "sync",
        module_root,
        f"s3://{bucket}/{args.module}",
        "--exclude",
        "*",
    ]
    if args.sync_results:
        sync_cmd += ["--include", "results/*"]
    if args.sync_plots:
        sync_cmd += ["--include", "plots/*"]
    # exclude hidden files (.dotfiles)
    sync_cmd += ["--exclude", "**/.*"]

    if args.delete_extra:
        sync_cmd += ["--delete"]
    if args.dryrun:
        sync_cmd += ["--dryrun"]
    if args.profile:
        sync_cmd += ["--profile", args.profile]

    sync_result = subprocess.run(sync_cmd)
    if sync_result.returncode:
        print(
            f"Error syncing to S3 bucket '{bucket}'.\n"
            "Check that the bucket name is correct and that"
            " you have the correct profile active (or use the --profile option)."
            " Be sure to have run `aws sso login` before running this script.\n",
            file=sys.stderr,
        )
        sys.exit(1)

    print("Sync complete.")
    print(f"Files are now available on S3 at: s3://{bucket}/{args.module}/")


if __name__ == "__main__":
    main()
