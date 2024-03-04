This file contains information about contributing to docs development.

`OpenScPCA-analysis` documentation is written in markdown with [`Material for MkDocs`](https://squidfunk.github.io/mkdocs-material/).
The [`mkdocs.yml`](../mkdocs.yml) is the documentation configuration file, and documentation pages are stored in this directory (`docs`).
To locally test documentation, you will need `python3` and dependencies listed in the [`requirements.txt`](./requirements.txt) file.
To install these dependencies, use your preferred approach:

```sh
# without conda
pip install -r requirements.txt

# with conda
conda install --file requirements.txt
```

To serve the documentation during local development, run `mkdocs serve` from the root directory of this repository.


## Docs organization

Documentation is written as a series of markdown files nested by topic in directories.

- All directories at the top-level inside `docs` will represent navbar sections.
  - Each top-level directory should contain an `index.md` file with an overall description of what that section contains.
  - The `index.md` should have an L1 header with the same title as the navbar section.
  - Markdown files in each directory represent sections shown along the left sidebar.
- Nested directories within each navbar section should be used to add an additional bold (gray) header on the left sidebar.
  - Markdown files within each nested directory will fall under this additional bold (gray) header.
- Any visual aids used in the docs should be placed in `docs/img`.

## Adding new documentation

As outlined in [this IA issue](https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/61), we have several sections of documentation.
The diagram below shows the directory structure for the docs, for now considering only the top three directory levels.
The directory names as written in the IA issue are shown in parentheses.

Each bold bullet point in the IA issue is a directory, and each plain-text bullet is a markdown file.
Consult the IA issue and this directory structure to determine where to place new markdown file.

If the IA indicates this file should be in an even more nested directory, don't worry about it for now; later, we will add additional nested directores as more documentation comes together.

Files and directories should be named with all lowercase letters and use `-`, not `_`, where needed.
If you are adding new images/visual aids for documentation, place them in the `img` directory (not shown in the diagram below).
**Be sure to consult the [style guide](./general-style-guide.md) when writing your documentation.**

```
├── welcome (Welcome to OpenScPCA)
│   ├── accessing-resources (Getting access to resources)
│   └── getting-started (Getting Started)
├── technical-setup (Setting up)
│   └── environment-setup (Getting access to resources)
├── communications-tools (Communicating within the Project)
│   ├── github-discussions (Github Discussions)
│   └── github-issues (Github Issues)
│       └── writing-issues (Writing good issues)
├── contributing-to-analyses (Contributing to an analysis)
│   ├── analysis-modules (Analysis modules)
│   ├── doing-analyses (Doing an analysis)
│   │   └── working-with-git (Working with git)
│   ├── planning-analyses (Planning your analysis)
│   └── pull-requests (Creating Pull Requests)
├── software-platforms (Software tools)
│   ├── aws
│   ├── docker
│   └── lsfr
└── troubleshooting-faq (Getting help & FAQ)
```
