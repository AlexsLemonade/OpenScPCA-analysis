# Running an analysis module

All modules should contain [clear documentation in the `README.md` file](documenting-analysis.md) about how to run them, including:

- Information about software and compute requirements
- What command(s) to issue to run the full module

!!! note
    If you are exploring modules and come across an analysis module without clear instructions for how to run it, please let us know using [GitHub Discussions](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/new?category=modify-an-existing-analysis).


In some cases, a module will only contain a single script (or notebook), and running that script alone will run the full analysis.

However, in most cases, modules will have several scripts (or notebooks) which need to be run in order, potentially each with specific input files.
If you are writing such a module, we strongly recommend you add a separate script (e.g., a shell script), to run all commands in order.

!!! tip
    For examples of a shell script that runs a module, please see those scripts in [our example R module](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/analyses/hello-R/run_hello-R.sh) or [our example Python module](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/analyses/hello-python/run_hello-python.sh).

You should save your script in the root of your analysis module and name it `run_{module-name}.sh` (this file name not strictly required, but it is highly recommended).

We also strongly recommend starting your script out with a few lines to make it easier to diagnose errors that might occur, and generally run more smoothly:

```bash
#!/bin/bash
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"
```

- `#!/bin/bash` is known as a ["shebang" line](https://linuxhandbook.com/shebang/).
It makes the script executable with a default bash interpreter.
- `set -euo pipefail` sets a few helpful script settings with the (you guessed it!) `set` command:
    - The `-e` flag will cause the script to fail immediately if any command it calls fails.
    - The `-u` flag will cause the script to fail immediately if it encounters an undefined variable.
    - The `-o pipefail` flag extends the `-e` flag to apply to commands that are part of pipelines.
- `cd "$(dirname "${BASH_SOURCE[0]}")"` automatically sets the working directory to the directory where the script is saved.
    - This means it will behave as though it is run from its own directory, no matter where you call the script from.
    - This is especially useful to make sure that all paths used within the script can be relative paths that will work no matter how the script is run.


After these lines, you can then proceed to add all of the commands and scripts needed to run the full module.
Be sure to add plenty of comments throughout the script to help orient other script users and contributors!
