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
        description="Sync daata from an analysis module to your S3 bucket.",
    )
    parser.add_argument(
        "--module",
        "-m",
        type=str,
        required=True,
        help="The name of the analysis module to sync data from.",
    )
    parser.add_argument(
        "--bucket",
        "-b",
        type=str,
        required=True,
        help="The name of the S3 bucket to sync data to.",
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
        "--no-delete",
        action="store_false",
        dest="delete_extra",
        help="Do not delete extra files in the S3 bucket.",
    )
    parser.add_argument(
        "--dryrun",
        action="store_true",
        help="Perform a dry run of the download: show what would be done but do not sync anything.",
    )
    parser.add_argument(
        "--profile",
        type=str,
        default=None,
        help="The AWS profile to use for the download. Uses the current default profile if undefined.",
    )

    args = parser.parse_args()

    ### Validate the arguments ###

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
        f"s3://{args.bucket}/{args.module}",
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
            f"Error syncing to S3 bucket '{args.bucket}'."
            " Check that the bucket exists and that you have permission to write to it.",
            file=sys.stderr,
        )
        print(sync_result.stderr, file=sys.stderr)
        sys.exit(1)

    print("Sync complete.")


if __name__ == "__main__":
    main()
