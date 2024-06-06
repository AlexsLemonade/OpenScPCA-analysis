# Structuring analysis scripts

While [notebooks](notebook-structure.md) are a great way to write analyses meant for interactive viewing and sharing, you may instead choose to use scripts (`*.R` or `*.py` for R and Python, respectively) to handle the "heavier lifting" aspects of your analysis.
Scripts are always a great place to add code or functions you tend to use frequently.

This page covers our recommendations for how you can structure your analysis scripts, which can store in the [`scripts` folder of your analysis module](./index.md#skeleton-analysis-module-contents).

To learn more about writing a wrapper script to run all scripts and/or notebooks in a module, please refer to [this documentation on running analysis modules](running-a-module.md).


Similar to our [recommendations for computational notebooks](notebook-structure.md), we recommend using a similar overall structure in your scripts to enhance sharing and reproducibility:

- **[Shebang](https://linuxhandbook.com/shebang/):** Add an opening line to set the script's interpreter.
- **Opening comments**: Add comments to the top of your script that briefly explain the script.
- **Load packages**: Load required packages.
- **Functions**: Define any custom functions used in the script.
All functions should include documentation comments.
- **Parse script arguments**: Parse any arguments to the script, such as file paths, here.
- **Main body**: This is where the main body of your script should appear.
- **Session info** (optional): Export out the versions of all packages used in the analysis.

Many considerations for writing scripts are the same as for writing notebooks, so we encourage you to read the documentation on structuring [R Markdown notebooks](notebook-structure.md#r-markdown-notebooks) and [Jupyter notebooks](notebook-structure.md#jupyter-notebooks).

Below, we provide a few additional considerations that are specific to scripts in general, and then specific considerations for R and Python.

### Considerations for any kind of script

- The first line of your R script should be a ["shebang,"](https://linuxhandbook.com/shebang/) which stands for "ha**sh** (`#`) **bang** (`!`)".
    - This makes the script executable by giving the exact path to the interpreter to use when running the script.
    - Below, we provide specific instructions for how to write this line in [R scripts](#r-shebang) and [Python scripts](#python-shebang).
- Your opening comments should be a brief, 1-3 sentences describing the purpose of the script.
    - These comments _are complementary to_ but do not replace [additional documentation in your module's `README.md`](documenting-analysis.md).
- Don't forget to document any custom functions that you write!
    - At a minimum, function documentation comments should describe the function's purpose, its arguments, and what it returns
- If you need to import a separate script with reusable functions, you should import it in the section where you load other packages.
- In most cases, your script should not hard code input and output file paths.
Instead, we recommend using an argument parser to read these paths in from the command line.
    - Below, we provide recommendations for argument parsing in [R scripts](#r-arguments) and [Python scripts](#python-arguments).
 - You may still want to export session information at the end of your script.
See our additional documentation on how to do this for [R scripts](../../ensuring-repro/managing-software/reporting-dependencies.md#in-r-scripts) and [Python scripts](../../ensuring-repro/managing-software/reporting-dependencies.md#in-python-scripts).

### Considerations for R scripts

#### Shebang { #r-shebang }

Your R scripts should always start with this line:

```r
#!/usr/bin/env Rscript
```

#### Parse input arguments { #r-arguments }

We recommend using the [`optparse` package](https://cran.r-project.org/package=optparse) to parse input arguments in R scripts.

For example, after loading the `optparse` package, you can define input options with the [`make_option()`](https://rdrr.io/cran/optparse/man/add_make_option.html) function, and then use the [`parse_args()`](https://rdrr.io/cran/optparse/man/parse_args.html) function to parse them into a usable list.

```r
# Set up arguments
option_list <- list(
  make_option(
    # Users can specify this argument with either `-i` or `--input_rds_file`
    opt_str = c("-i", "--input_rds_file"),
    # The provided value will be stored as a character variable called `input_rds_file`
    type = "character",
    dest = "input_rds_file",
    # Description of the argument
    help = "path to input RDS file for analysis"
  ),
  ## define any additional options here ##
)

# Parse arguments into a usable list
opts <- parse_args(OptionParser(option_list = option_list))
```

In the above example, you can access the parsed argument as `opts$input_rds_file`.

#### Additional R tips

- You can add more structure to your script using [sectioning comments](https://r4ds.hadley.nz/workflow-style#sectioning-comments), which help you navigate the code in RStudio much like section headers in a notebook.
- Please use [docstrings](https://www.geeksforgeeks.org/python-docstrings/) when adding documentation to any custom functions you write.


### Considerations for Python scripts

#### Shebang { #python-shebang }

Your Python scripts should always start with this line:

```python
#!/usr/bin/env python3
```

#### Parse input arguments { #python-arguments }

We recommend using the [`argparse` package](https://docs.python.org/3/library/argparse.html) to parse input arguments in Python scripts.

For example, after loading the `argparse` package, you can define a parser object with [`argparse.ArgumentParser()`](https://docs.python.org/3/library/argparse.html#creating-a-parser), and then use the [`.add_argument()` method](https://docs.python.org/3/library/argparse.html#adding-arguments) to add each argument to the parser object.
Finally, use the [`.parse_args()` method](https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser.parse_args) to create an object containing the argument values.


```python
# define a parser object with helpful information for script users
parser = argparse.ArgumentParser(
    description = "Performs task XYZ as part of analysis"
)

# define a command line argument
parser.add_argument(
    # Users can specify this argument with either `-i` or `--input_anndata_file`
    "-i",
    "--input_anndata_file",
    # The provided value will be stored as a string called `input_anndata_file`
    type = "character",
    dest = "input_adddata_file",
    # Description of the argument
    help = "path to input AnnData file for analysis"
)
## define any additional arguments here ##

# parse the arguments into an object
args = parser.parse_args()
```

In the above example, you can access the parsed argument as `args.input_anndata_file`.



#### Additional Python tips

- Include [docstrings](https://www.geeksforgeeks.org/python-docstrings/) for all custom functions you write.
    - You may wish to enclose your opening script comments inside a docstring.
- Following Python convention, you may also wish to define [a `main()` function](https://realpython.com/python-main-function/) to handle most of the script work.
    - In this case, your [argument parsing code](#python-arguments) should be placed at the top of the `main()` function.
