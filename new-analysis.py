#!/usr/env/bin python3

""" Script to create a new analysis module for OpenScPCA """

import argparse
import os
import pathlib
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(
        prog="new-analysis", description="Create a new analysis module for OpenScPCA."
    )
    parser.add_argument(
        "name",
        help=(
            "Name of the analysis module."
            " This will be used for the directory name and must not already exist within the `analyses` directory."
        ),
    )
    parser.add_argument(
        "--use-conda",
        action="store_true",
        default=False,
        help="Set up a new conda environment file for the module.",
    )
    parser.add_argument(
        "--use-renv",
        action="store_true",
        default=False,
        help="Initialize a new renv environment for the module.",
    )
    parser.add_argument(
        "--condaenv-name",
        help=(
            "Name of the conda environment to create for the module."
            " The default name is `openscpca_<name>` where <name> is the module name."
        ),
    )

    args = parser.parse_args()

    # get the paths relative to this script file
    base_dir = pathlib.Path(__file__).parent
    template_dir = base_dir / "templates" / "analysis-module"
    module_dir = base_dir / "analyses" / args.name

    # exit if the module directory already exists
    if module_dir.exists():
        print(
            f"Analysis module `{args.name}` already exists at `{module_dir}`.",
            "Exiting.",
            sep=os.linesep,
        )
        return 1

    # if use_conda is requested, check that conda is available
    if args.use_conda and not shutil.which("conda"):
        print(
            "Setup with conda was requested, but conda is not available on the system.",
            "Please install conda (Miniconda) and try again.",
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
    try:
        shutil.copytree(template_dir, module_dir)
    except FileNotFoundError:
        # just in case the template directory is missing
        if not template_dir.exists():
            print("Expected template directory does not exist at `{template_dir}`.")
            print("Exiting.")
            return 1
        else:
            raise

    if args.use_conda:
        # add an environment.yml file copied from base but with a new name based on the module name
        if args.condaenv_name:
            env_name = args.condaenv_name
        else:
            env_name = f"openscpca_{args.name}"

        with open(base_dir / "environment.yml", "r") as f:
            lines = f.readlines()
        with open(module_dir / "environment.yml", "w") as f:
            for line in lines:
                if line.startswith("name:"):
                    f.write(f"name: {env_name}\n")
                else:
                    f.write(line)

        # should we run conda env create -f environment.yml here?
        # if we do, we need to check that the environment name does not already exist

        # subprocess.run(
        #     ["conda", "env", "create", "-f", "environment.yml"],
        #     cwd=module_dir,
        # )

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
