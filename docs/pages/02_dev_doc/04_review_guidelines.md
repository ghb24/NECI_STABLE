---
title: Review guidelines
---

## Review guidelines

Any new code which is to be added to the master branch of NECI has to
undergo code review to check if it does not introduce new bugs into
existing code and if the code is written in an understandable and
maintanable way, and following the code conventions introduced in
section
<a href="#sec:code_conv" data-reference-type="ref" data-reference="sec:code_conv">1.1</a>.

The process of adding new code to the program contains these steps

1.  Create a new git branch and add your code there. Do not forget to
    add tests, so the functionality can be verified.

2.  Push the branch to the bitbucket repository

3.  Create a pull request to master/devel (depending on where you want
    the code merged), selecting a set of reviewers.

4.  The pull request triggers a pipeline that tries to compile the
    program and run a set of tests with the new code. If any of those
    fail, you will get a notification and you can check the reason
    therefore.

5.  The reviewers will now check the code and can comment on it. Make
    sure to address these comments.

6.  Once you got the approval of at least one reviewer, you can merge
    the code. It is now included in the program.

If you have been selected as a reviewer and decide to do the review,
check the pull request. It will contain information on the ran tests and
a list of all changes to the code. Go through the new code and check if
it is written in a clean and well-commented way, in accordance with the
code conventions. Does it have tests? You do not have to verify that the
new code is bug-free nor do you need to debug it, that is not within the
scope of the review. If you find something that could be improved, make
a comment on that or create a task, this helps the author to increase
the quality of the code.
