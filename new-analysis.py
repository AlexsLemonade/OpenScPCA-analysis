#!/usr/env/bin python3

import argparse
import os
import pathlib
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(
        prog="new-analysis", description="Create a new analysis module for OpenScPCA"
    )
    parser.add_argument(
        "name",
        help=(
            "Name of the analysis module."
            " This will be used for the directory name and must not already exists within the `analyses` directory."
        ),
    )
    parser.add_argument(
        "--use-conda",
        action="store_true",
        default=False,
        help="Set up a new conda environment file for the module",
    )
    parser.add_argument(
        "--use-renv",
        action="store_true",
        default=False,
        help="Initialize a new renv environment for the module",
    )

    args = parser.parse_args()

    # get the paths relative to this script file
    base_dir = pathlib.Path(__file__).parent
    template_dir = base_dir / "analyses" / "_template"
    module_dir = base_dir / "analyses" / args.name

    # exit if the directory already exists
    if module_dir.exists():
        print(
            f"Analysis module `{args.name}` already exists at `{module_dir}`.",
            "Exiting.",
            sep=os.linesep,
        )
        return 1

    # if use_renv is requested, check that R is available
    if args.use_renv and not shutil.which("Rscript"):
        print(
            "Setup with renv was requested, but Rscript is not available on the system.",
            "Please install R and try again.",
            sep=os.linesep,
        )
        return 1

    # create the new module directory from the template
    shutil.copytree(template_dir, module_dir)

    if args.use_conda:
        # add an environment.yml file copied from base but with the new name
        with open(base_dir / "environment.yml", "r") as f:
            lines = f.readlines()
        with open(module_dir / "environment.yml", "w") as f:
            for line in lines:
                if line.startswith("name:"):
                    f.write(f"name: {args.name}\n")
                else:
                    f.write(line)

        # should we run conda env create -f environment.yml here?

    if args.use_renv:
        # initialize a new renv environment
        renv_script = """
        if (!requireNamespace("renv", quietly = TRUE))
            install.packages("renv")
        renv::scaffold()
        """

        subprocess.run(
            ["Rscript", "-e", renv_script],
            cwd=module_dir,
        )


if __name__ == "__main__":
    main()
