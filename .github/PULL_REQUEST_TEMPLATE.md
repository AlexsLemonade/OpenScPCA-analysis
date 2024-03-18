## For PRs updating docs, please use the appropriate template:

Use the 'Preview' view to click the link below:

<a href="?expand=1&template=docs-update-pr.md"> Update docs PR template </a>

### For any other types of PRs, delete all of the above text and this line and use the template below.


_Hi there, thanks for your contribution!_
_Please fill out this template to help us review your code._

_Remember, the more context you provide, the faster the review will go!_
_For more information about filling out this template, please see the [OpenScPCA documentation on filing pull requests](STUB_LINK docs on PRs)._

_Before filing the pull request, you can also feel free to delete the italicized instruction lines in this template._



### Purpose/implementation Section

_In this section, tell reviewers about the scope and purpose of the code in the pull request._


#### What is the goal of this pull request?

_Include the scientific or analytical goals, if applicable._



#### Briefly describe the general approach you took to achieve this goal.




#### Please link to the GitHub issue that this pull request addresses.


#### If known, do you anticipate filing additional pull requests to complete this analysis module?



### Results

_Delete this section if no results are associated with your pull request._

_In this section, tell reviewers about what kinds of results (if any) your code produces._
#### What is the name of your results bucket on S3?

_Results should be [uploaded to your bucket](STUB_LINK docs on uploading results to S3) so they are available during review._


#### What types of results does your code produce (e.g., table, figure)?


#### What is your summary of the results?




### Provide directions for reviewers

_In this section, tell reviewers what kind of feedback you are looking for._
_This information will help guide their review._

#### Are there particularly areas you'd like reviewers to have a close look at?



#### Is there anything that you want to discuss further?






### Author checklists

_Check all those that apply._
_Note that you may find it easier to check off these items after the pull request is actually filed._


#### Analysis module and review

- [ ] This analysis module [uses the analysis template and has the expected directory structure](STUB_LINK for analysis structure doc).
- [ ] The analysis module `README.md` has been updated to reflect code changes in this pull request.
- [ ] The analytical code is documented and contains comments.
- [ ] Any results this code produces have been added to S3 for review.

#### Reproducibility checklist

- [ ] Code in this pull request has been added to the GitHub Action workflow that runs this module.
- [ ] The dependencies required to run the code in this pull request have been added to the analysis module `Dockerfile`.
- [ ] If applicable, the dependencies required to run the code in this pull request have been added to the analysis module conda `environment.yml` file.
- [ ] If applicable, R package dependencies required to run the code in this pull request have been added to the analysis module `renv.lock` file.
