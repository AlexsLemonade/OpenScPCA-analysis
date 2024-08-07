This file contains information about contributing to docs development.

`OpenScPCA-analysis` documentation is written in markdown with [`Material for MkDocs`](https://squidfunk.github.io/mkdocs-material/).
The `../mkdocs.yml` file is the documentation configuration file, and documentation pages are stored in this directory (`docs`).
To locally test documentation, you will need `python3` and dependencies listed in the [`requirements.txt`](./requirements.txt) file.
To install these dependencies, use your preferred approach:

```sh
# without conda
pip install -r requirements.txt

# with conda
conda install --file requirements.txt
```

To serve the documentation during local development, run `mkdocs serve` from the root directory of this repository.

**Table of contents**
- [Docs organization](#docs-organization)
- [Adding new documentation](#adding-new-documentation)
- [Linking to other pages](#linking-to-other-pages)
- [Adding images](#adding-images)
- [Phrasing and other syntax](#phrasing-and-other-syntax)


## Docs organization

Documentation is written as a series of markdown files nested by topic in directories.

- All directories at the top-level inside `docs` will represent `navbar` sections.
  - Each top-level directory should contain an `index.md` file with an overview of what that section contains.
  - The `index.md` should have an L1 header with the same title as the `navbar` section.
  - Markdown files in each directory represent sections shown along the left sidebar.
- Nested directories within each `navbar` section should be used to add an sections along the left sidebar.
  - There should only be _one level_ of nested directories (i.e., no directories inside these ones)
  - Nested directories may also have an `index.md` file with an L1 header and an overall description of what that section contains.
- Any visual aids used in the docs should be placed in `docs/img` (see the [adding images section](#adding-images)).

## Adding new documentation

As outlined in [this IA issue](https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/61), we have several sections of documentation.
The diagram below shows the directory structure for the docs.
The directory names as written in the IA issue are shown in parentheses.

Each bold bullet point in the IA issue is a directory, and each plain-text bullet is a markdown file.
Consult the IA issue and this directory structure to determine where to place new markdown file.

If the IA indicates this file should be in an even more nested directory, don't worry about it for now; later, we will add additional nested directories as more documentation comes together.

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
│   └── github-issues (Github Issues)
│       └── writing-issues (Writing good issues)
├── contributing-to-analyses (Contributing to an analysis)
│   ├── analysis-modules (Analysis modules)
│   ├── doing-analyses (Doing an analysis)
│   │   └── working-with-git (Working with git)
│   ├── planning-analyses (Planning your analysis)
│   ├── pull-requests (Creating Pull Requests)
│   │   └── filing-pull-requests (Filing a PR)
│   │       ├── before-filing (Before filing a PR)
│   │       └── creating-pull-requests (Creating the PR)
│   ├── pr-review-process (PR Review Process)
│   └── merging-prs (Merging a PR)
├── software-platforms (Software tools)
│   ├── aws
│   ├── docker
│   └── lsfr
└── troubleshooting-faq (Getting help & FAQ)
```

## Linking to other pages

It is helpful to include relative link to other docs pages when discussing additional concepts beyond the scope of the given docs page.
In many circumstances during docs development, those pages will not yet have been written.
In those circumstances, use the phrase `#STUB_LINK` where you would otherwise add a relative link.
For example:

```
This sentence should link to the [pull request documentation](#STUB_LINK).
```


## Adding images

All image files should be placed in `docs/img`.
To include an image in a docs file, use the following syntax to center the image on the page:

```
    <figure markdown="span">
      ![alt text](relative link to image file){width="number of pixels"}
    </figure>
```

When adding the image file itself, try to make sure that its width is roughly 2X the number of pixels you specify in the docs.
For example, you might follow this procedure:

- Create an image file (don't worry about the size yet)
- Add the image to your docs file and determine a good pixel width for the embedded figure with some trial/error
- If needed, resize the image file itself to be 2X the pixel width you specified in the embedded figure

We have found that an embedded figure width of ~400-600 pixels seems to work well.
In this case, the associated image file would be ~800-1200 pixels wide.


## Phrasing and other syntax

Only these items should be in backticks:
- Actual code
- Names of packages, e.g. `renv`
- SemVar software versions
- References to a specific repository, e.g. `OpenScPCA-analysis`
- Buttons in GitKraken, e.g.. `Stage all files`

For consistency, please use the following:

- Refer to `github.com` as "GitHub" (but without quotes) except when providing specific links.
- Refer to git as "Git" (without quotes), not "git"
  - If you are actually referring to a git command, you should use `git <command>`
- Capitalize Python (not python)
- Capitalize Docker (not docker)
