---
name: Create a release
about: Use this to track creating a semantically-versioned release
title: 'Create a release:'
labels: release
---

If issues must be resolved before creating a release, mark them as blockers in ZenHub.

**Pre-release checklist:**

- [ ] Any mentions of the release tag in the repository have been updated
- [ ] All scheduled workflows are passing
- [ ] All blocking issues have been resolved
- [ ] Write release notes, which should include the following:
      - Which modules have been added, if any?
      - Which modules have been deprecated, if any?
      - What has changed in the repo documentation, if anything?
      - Are there any _major_ shifts in project dependencies?
      For example, is there a package that was used throughout Docker environments that has been replaced?
      - Have there been any changes in repo file organization?
