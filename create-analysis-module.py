#!/usr/bin/env python3

"""Script to create a new analysis module for OpenScPCA"""

import argparse
import pathlib
import re
import shutil
import subprocess
import sys


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a new analysis module for OpenScPCA.",
        usage="%(prog)s MODULE_NAME [options]",
    )
    parser.add_argument(
        "name",
        metavar="MODULE_NAME",
        help=(
            "Name of the analysis module to create."
            " This will be used for the directory name and must not already exist within the `analyses` directory."
        ),
    )
    parser.add_argument(
        "--use-renv",
        action="store_true",
        default=False,
        help="Initialize a new renv environment for the module.",
    )
    parser.add_argument(
        "--use-conda",
        action="store_true",
        default=False,
        help="Initialize a new conda environment for the module.",
    )
    parser.add_argument(
        "--conda-file-only",
        action="store_true",
        help=(
            "Only make the conda `environment.yml` file."
            " Creating a new conda environment is skipped."
        ),
    )

    args = parser.parse_args()

    # get the paths relative to this script file
    base_dir = pathlib.Path(__file__).parent
    template_dir = base_dir / "templates" / "analysis-module"
    module_dir = base_dir / "analyses" / args.name

    # fail if the module name is not a simple directory name
    if not re.search(r"^[A-Za-z0-9_\-]+$", args.name):
        sys.exit(
            "Module name should not contain spaces or special characters.\nExiting."
        )

    # exit if the module directory already exists
    if module_dir.exists():
        sys.exit(
            f"Analysis module `{args.name}` already exists at `{module_dir}`."
            "\nExiting."
        )

    # set the conda environment name
    env_name = f"openscpca-{args.name}"

    # if use_conda is requested, check that conda is available
    if args.use_conda and not args.conda_file_only:
        if not shutil.which("conda"):
            sys.exit(
                "Setup with conda was requested, but conda is not available on the system."
                "\nPlease install conda (Miniconda) and try again."
            )
        # check for existing environment name
        existing_envs = subprocess.run(
            ["conda", "env", "list"],
            capture_output=True,
            text=True,
        ).stdout
        existing_envs = [env.split()[0] for env in existing_envs.splitlines() if env]
        if env_name in existing_envs:
            sys.exit(
                f"\nConda environment `{env_name}` already exists."
                "\nYou may use the `--conda_file_only` option to create the environment file only and manually create your environment with:"
                "\n`conda env create --f analyses/<module>/environment.yml --name <env_name>`"
                "\n\nAlternatively, you can remove the existing environment with:"
                f"\n`conda env remove --name {env_name}`"
                "\nor rename the existing environment:"
                f"\n\n`conda rename --name {env_name} <new_name>`"
                "\n\nExiting."
            )

    # if use_renv is requested, check that R is available
    if args.use_renv and not shutil.which("Rscript"):
        sys.exit(
            "Setup with renv was requested, but Rscript is not available on the system."
            "\nPlease install R and try again."
        )

    # create the new module directory from the template
    try:
        shutil.copytree(template_dir, module_dir)
        print(f"\nCreated new analysis module in `{module_dir}`.")
    except FileNotFoundError:
        # just in case the template directory is missing
        if not template_dir.exists():
            sys.exit(
                "Expected template directory does not exist at `{template_dir}`."
                "\nExiting."
            )
        else:
            raise

    if args.use_conda or args.conda_file_only:
        # add an environment.yml file copied from base but with a new name based on the module name
        with open(base_dir / "environment.yml", "r") as f:
            lines = f.readlines()
        with open(module_dir / "environment.yml", "w") as f:
            for line in lines:
                if line.startswith("name:"):
                    f.write(f"name: {env_name}\n")
                else:
                    f.write(line)
        print(
            f"\nCreated a conda environment file `environment.yml` in `{module_dir}`."
        )

        # create the conda environment
        if not args.conda_file_only:
            subprocess.run(
                ["conda", "env", "create", "-f", "environment.yml"],
                cwd=module_dir,
            )
            print(f"Created conda environment `{env_name}`.")

            # alternatively, we could create the environment in a module subdirectory with the `--prefix` option
            # subprocess.run(
            #     ["conda", "env", "create", "-f", "environment.yml", "--prefix", "env"],
            #     cwd=module_dir,
            # )
            # with open(module_dir / "env" / ".gitignore", "w") as f:
            #     f.writelines("# ignore all environment files\n*\n")

    if args.use_renv:
        # initialize a new renv environment
        renv_script = """
            if (!requireNamespace("renv", quietly = TRUE))
                install.packages("renv")
            renv::scaffold(
                repos = list(CRAN = "https://p3m.dev/cran/latest"),
                settings = list(
                    ppm.enabled = TRUE,
                    r.version = "4.3.3",
                    bioconductor.version = "3.18"
                )
            )
        """

        subprocess.run(
            ["Rscript", "-e", renv_script],
            cwd=module_dir,
        )
        # make the components directory and add a dependencies.R file
        component_dir = module_dir / "components"
        component_dir.mkdir(exist_ok=True)
        (component_dir / "dependencies.R").write_text(
            "# R dependencies not captured by `renv`\n" '# library("missing_package")\n'
        )

        print(f"\nInitialized new renv environment in `{module_dir}`.")


if __name__ == "__main__":
    main()
