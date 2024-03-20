# Analysis modules

OpenScPCA organizes individual analyses into [analysis modules](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses).
Each analysis module is a folder with files containing all code, computing environment specifications, and documentation needed to run and interpret an analysis.
For example, an analysis to perform cell type annotation on Ewing sarcoma samples would be single analysis module named, and it might be named `celltype-ewings-sarcoma`.

This section explains the structure of analysis modules.

## Skeleton analysis module contents

You can create a starting point for your analysis module with the script `create-analysis-module.py`.
This script will create a skeleton analysis module with the following file structure:

<figure markdown="span">
    ![Baseline file structure of an analysis module.](STUB_LINK to image created in #134){width="200"}
</figure>

Please refer to [the documentation on creating analysis modules](STUB_LINK for creating modules) when you are ready to make your first analysis folder and begin contributing to OpenScPCA!

These are the main files and folders you will interact with when writing your analysis:

- **`scripts`**
    - You can save any scripts (e.g., `.R`, `.py`, or `.sh`) that you write for your analysis module in this folder.
    - If you choose, you can also save any notebooks (e.g., R Markdown or Jupyter) files in this folder too.
    Or, depending on the goals of your analysis module, you may prefer to save notebooks in the top-level of the analysis module folder.
- **`results`**
    - Any result files (e.g., TSV files) that your code produces should be saved to this `results` folder.
    - Git will ignore the contents of this folder, _except_ for its `README.md` file, which you can use to document the results files themselves.
    This means that only its `README.md` file will be present in the remote repository.
- **`plots`**
    - Any plots that your code produces should be saved to this `plots` folder.
- **`scratch`**
    - You can optionally use this folder to store _intermediate_ files that your code produces but are not meant to live in `results` or `plots`.
    - Git will ignore the contents of this folder, so contents of this folder will only be stored locally and not in the remote repository.
- **`README.md`**
    - Use this [markdown file](STUB_LINK docs on markdown) to document your analysis module.
  Your `README.md` file should have enough information for other contributors or repository users to learn the following:
        - The scientific goals of the module
        - The input data the module uses and its [computational resource requirements](../starting-your-analysis/determining-compute-requirements.md)
        - How to run the module
    - Please see the documentation on [documenting your analysis module](STUB_LINK module documenting notebooks) for more information about adding to this `README.md` file.


There are also some additional files in the skeleton that are useful to be aware of:


- **`Dockerfile`**
    - This is the analysis module's [Dockerfile](https://docs.docker.com/reference/dockerfile/) and contains the commands that Docker uses to build the module's Docker image.
    - For more information on how OpenScPCA uses `Docker` images, [please see our `Docker` documentation](STUB_LINK for docker docs).
- Hidden files **`.gitignore`** and **`.dockerignore`**
    - We have set up these files to tell Git and Docker, respectively, to ignore certain files that do not belong in version control or in the module's Docker image.
    - These files will likely be automatically hidden from you, and you don't really have to worry about it.
    Just be aware that they are there and working behind the scenes to help manage the module!



## Additional files you will add to your module

While you write your analysis, you may add other files too:

- Scripts and analysis notebooks, e.g., R Markdown files or Jupyter notebooks
    - We recommend saving scripts in the `scripts` folder, as described above.
    - You are also welcome to notebook files in `scripts`, or in the top-level of your analysis module folder.
    Feel free to choose what location makes the most sense for your analysis, as long as it is all documented in the module's `README.md` file!
        - Please see the documentation on [structuring your analysis notebooks](STUB_LINK structure notebooks) for more information about how to write your analysis notebooks.

    !!! note "Naming your files"
        If your module has multiple scripts or notebooks, we recommend naming them in the order they should be run.

        For example, you might have these files in your module (the script names are conceptual!):

        - `scripts/01_script-to-prepare-data.R`
        - `scripts/02_script-to-analyze-data.R`
        - `03_notebook-to-visualize-data.Rmd`

- A shell script to run all code in the module
    - If your module has multiple scripts or notebooks, we recommend adding a script (for example a `shell` script) to the top-level of your analysis module folder that will run all scripts in order.
    - You can name this file, for example, `run_<module name>.sh` and document how to run in the module's `README.md`
- Additional environment files
    - When you [create a module](STUB_LINK for creating a module), you can choose to include files that manage the module's software environment in the analysis module skeleton.
    - In this case, your module may also contain contain R-specific (e.g., `renv.lock`) and/or Python-specific (e.g., `environment.yml`) files or folders.




## Example analysis modules

To help you get started, we've created two analysis modules that you can use as references:

- [`hello-R`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/hello-R) is an example R-based analysis module
- [`hello-python`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/hello-python) is an example Python-based analysis module
