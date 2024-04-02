# Structuring analysis scripts

While [notebooks](notebook-structure.md) are a great way to write analyses meant for interactive viewing and sharing, you may instead choose to use scripts (`*.R` or `*.py` for R and Python, respectively) to handle the "heavier lifting" aspects of your analysis.
Scripts are always a great place to add code or functions you tend to use frequently.

This page covers our recommendations for how you can structure your analysis scripts.

To learn more about writing a wrapper script to run all scripts and/or notebooks in a module, please [see this documentation on running analysis modules](running-a-module.md).

## Saving your scripts

You can store these in the [`scripts` folder of your analysis module](./index.md#skeleton-analysis-module-contents).
If you find yourself writing utility functions or scripts that do not necessarily represent the main analysis code itself, you are welcome to add a `utils` folder in your module to store these scripts.

As you develop your analysis, you may find yourself adding more than one script.
Often, these scripts are meant to be run in a particular order.
In addition to [documenting](documenting-analysis.md) the relationship of your scripts to one another, we recommend _naming_ your scripts with the overall module organization in mind.
Specifically, naming your scripts based on the order they are meant to be run can be helpful for you and other users, e.g. for an R-based module:

  - `scripts/01_first-module-script.R`
  - `scripts/02_second-module-script.R`
  - `scripts/03_third-module-script.R`
  - etc.

!!! tip "Naming files"
    Check out [these slides from Jenny Bryan](https://speakerdeck.com/jennybc/how-to-name-files) to learn more about how we in the Data Lab think about naming our files!


## Organizing your scripts

Similar to our [recommendations for computational notebooks](notebook-structure.md), we recommend using a similar overall structure in your scripts to enhance sharing and reproducibility:

<!-- Do we want to recommend a shebang? -->
- **Opening comments**: Add comments to the top of your script that briefly explain the script.
- **Load packages**: Load required packages.
- **Functions**: Define any custom functions used in the script.
- **File paths**: Define any input and output file paths.
- **Analysis**: This is where the bulk of the analysis should appear.
- **Session info**: Export out the versions of all packages used in the analysis.

Most considerations for writing scripts are the same as for writing notebooks, so we encourage you to read the documentation on structuring [R Markdown notebooks](notebook-structure.md#r-markdown-notebooks) and [Jupyter notebooks](notebook-structure/#jupyter-notebooks).

Below, we provide a few additional considerations that are specific to scripts in general, and then specific considerations for R and Python.

### Considerations for any kind of script

- Your opening comments should be a brief, 1-3 sentences describing the purpose of the script.
These comments _are complementary to_ but do not replace [additional documentation in your module's `README.md`](documenting-analysis.md).
- If you need to import a script with reusable functions, you should import it in the section where you load other packages.
- You will still need to export session information at the end of your script.
See our additional documentation on how to do this for [R scripts](../determining-requirements/determining-software-requirements.md#in-a-script) and [Python scripts](../determining-requirements/determining-software-requirements.md#in-python-scripts).


### Considerations for R scripts

- You can add more structure to your script using [sectioning comments](https://r4ds.hadley.nz/workflow-style#sectioning-comments), which help you navigate the code in RStudio much like section headers in a notebook.


### Considerations for Python scripts

- You may wish to enclose your opening comments inside a [triple-quote block](https://www.geeksforgeeks.org/triple-quotes-in-python/).
