# Code review and merging your pull request

Once you have submitted your [pull request (PR)](../creating-pull-requests/index.md), the Data Lab team will assign a reviewer to your pull request.
You can expect that each PR will undergo several rounds of code review and response to review.

Your reviewer will evaluate the code and/or documentation in your PR for the following considerations:

- **Context**
    - Has enough information been provided for the reviewer to fully understand the scope and context of what they are reviewing?
    If not, your reviewer's initial code review will likely request additional information they need to be able to perform review.
- **Clarity and Correctness**
    - Is the code readable, reasonably efficient, and well-commented?
    - Does the code appear to work correctly?
- **Documentation**
    - Are the code and results (if applicable) clearly documented?
    - Are the steps to set up the code environment and run the code clearly documented?
- **Reproducibility**
    - Can the code be re-run to successfully generate the same results?
- **Code checks**
    - We have set up several automated checks using GitHub actions for code quality control.
    If any of these automated checks fail, your reviewer may include feedback about what changes you may need to make for code checks to pass.


## Overview of the review process

Below is a general overview of what you can expect during the code review process:

1. In most cases, the first time your PR is reviewed by a Data Lab member, they will provide comments and suggestions. <!-- Feel free to browse some [example reviewer comments](STUB_LINK to example review comments). -->

1. You will then need to respond to those comments. Follow [these guidelines when addressing review comments](./respond-to-review.md).

1. Once you have fully addressed comments from the reviewer, then you will be able to re-request a review.

The above steps will repeat until the reviewer feels that all comments have been addressed adequately.
At that time, the PR will be approved.
Bear in mind that it is _normal and expected_ for a pull request to undergo multiple rounds of review before it is approved.

## Approved pull requests

Once approved, the PR will be merged into the main branch of `AlexsLemonade/OpenScPCA-analysis` by a Data Lab team member.

!!! tip

    After your PR has been merged, be sure to [sync your fork with `AlexsLemonade/OpenScPCA-analysis`](../working-with-git/staying-in-sync-with-upstream.md) to ensure you are working with the most up to date code.
