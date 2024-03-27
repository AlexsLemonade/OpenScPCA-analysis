#!/usr/bin/env python3

"""Script to create a new analysis module for OpenScPCA"""

import argparse
import pathlib
import re
import shutil
import subprocess
import sys
from typing import Union


def copy_file_with_tag_replacement(
    src: Union[pathlib.Path, str],
    dest: Union[pathlib.Path, str],
    tag: str,
    replacement: str,
) -> None:
    """
    Copy a file from src to dest, replacing all occurrences of tag with replacement.

    The tag in the file must be enclosed in `{{` and `}}` and may contain whitespace.
    """
    tag = tag.strip(
        " {}"
    )  # make sure the incoming tag does not have enclosing {} or whitespace
    tag_re = re.compile(r"{{\s*" + tag + r"\s*}}")
    with open(src, "r") as f:
        content = f.readlines()
    # replace the module tag with the module name
    output = (tag_re.sub(replacement, line) for line in content)
    with open(dest, "w") as f:
        f.writelines(output)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a new analysis module for OpenScPCA.",
        usage="./%(prog)s MODULE_NAME [options]",
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
        "-R",
        "--use-r",
        "--use-R",
        dest="use_r",
        action="store_true",
        default=False,
        help="Set up for R analysis, including an R Markdown notebook template.",
    )
    parser.add_argument(
        "--use-renv",
        action="store_true",
        default=False,
        help="Initialize a new renv environment for the module. Implies `--use-r`.",
    )
    parser.add_argument(
        "--use-conda",
        action="store_true",
        default=False,
        help="Initialize a new conda environment for the module.",
    )
    parser.add_argument(
        "--use-python",
        action="store_true",
        default=False,
        help="Set up for Python analysis, including a Python script template. Implies `--use-conda`.",
    )
    parser.add_argument(
        "--use-jupyter",
        action="store_true",
        default=False,
        help="Set up for analysis with Jupyter, including a Jupyter Notebook template. Implies `--use-conda`.",
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

    if args.use_renv:
        args.use_r = True

    if args.use_jupyter or args.use_python:
        args.use_conda = True

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

    # set up final status messages
    final_messages = [f"Creating analysis module `{args.name}` complete:\n"]

    # create the new module directory from the template
    try:
        shutil.copytree(template_dir, module_dir)
        final_messages.append(
            f"- Created new analysis module directory: `{module_dir}`."
        )
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
        # add a template environment.yml file
        if args.use_jupyter:
            env_template = base_dir / "templates" / "jupyter" / "environment.yml"
        else:
            env_template = base_dir / "templates" / "python" / "environment.yml"
        module_env = module_dir / "environment.yml"
        copy_file_with_tag_replacement(
            src=env_template,
            dest=module_env,
            tag="openscpca_module",
            replacement=args.name,
        )
        final_messages.append(f"- Created conda environment file: `{module_env}`.")

        # create the conda environment
        if not args.conda_file_only:
            subprocess.run(
                ["conda", "env", "create", "-f", "environment.yml"],
                cwd=module_dir,
            )
            final_messages.append(f"- Created conda environment: `{env_name}`.")

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

        final_messages.append(f"- Initialized new renv environment in `{module_dir}`.")

    # Add template files

    if args.use_r:
        template_rmd = base_dir / "templates" / "rmarkdown" / "notebook-template.Rmd"
        module_rmd = module_dir / "notebook-template.Rmd"
        copy_file_with_tag_replacement(
            src=template_rmd,
            dest=module_dir / "notebook-template.Rmd",
            tag="openscpca_module",
            replacement=args.name,
        )

        final_messages.append(f"- Added R Markdown notebook template: `{module_rmd}`.")

    if args.use_python:
        template_py = base_dir / "templates" / "python" / "script-template.py"
        module_py = module_dir / "script-template.py"
        copy_file_with_tag_replacement(
            src=template_py,
            dest=module_py,
            tag="openscpca_module",
            replacement=args.name,
        )

    if args.use_jupyter:
        template_ipynb = base_dir / "templates" / "jupyter" / "notebook-template.ipynb"
        module_ipynb = module_dir / "notebook-template.ipynb"
        copy_file_with_tag_replacement(
            src=template_ipynb,
            dest=module_ipynb,
            tag="openscpca_module",
            replacement=args.name,
        )

        final_messages.append(f"- Added Jupyter Notebook template: `{module_ipynb}`.")

    # print final status messages
    print()  # add a newline before the final messages
    print("\n".join(final_messages))


if __name__ == "__main__":
    main()
