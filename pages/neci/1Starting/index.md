title: Installation
---
# Using NECI

[TOC]
Git overview

It is essential if you plan to do developmental work to get familiar
with the source-code management software ‘git’. The code will get
unusable exponentially quickly if all development and new ideas are
hacked into the master branch of the code. The nature of research is
that most things probably won’t work, but you want to implement them and
test relatively quickly, without requiring a standard of code that will
remain usable in perpetuity. To avoid an inexorable increase in code
‘clutter’, it is essential to work in ‘branches’ off the main code. For
a more detailed introduction to the git package, see
[git-scm.com/book/en/v2/getting-started-git-basics](git-scm.com/book/en/v2/getting-started-git-basics).
In short, the workflow should be:

1.  Branch off a clean master version to implement something

2.  Test and develop in the branch

3.  Regularly merge the new code from the master branch into your
    personal development branch

4.  Once satisfied with the development, and that it is an improvement
    in scope or efficiency of the existing code, ensure it is tidy,
    commented, documented, as bug-free as possible, and tests added to
    the test suite for it. This may involve reimplementing it from a
    clean version of master if it can be done more efficiently

5.  Merge code back into master branch

A few potentially useful git commands in roughly the workflow described
above:

-   **git branch**  
    See what branch I am on. -a flag for all (inc. remote) branches.

-   **git pull origin master**  
    Update the master branch into the current local repository

-   **git checkout -b newbranchname**  
    Fork off current branch to a new branch called ‘newbranchname’

-   **git commit -a -m ‘Commit message’**  
    Commit a set of changes for the current branch to your local
    repository.

-   **git push origin branchname**  
    Push your current local branch called branchname to a new remote
    branch of the same name to allow access to others and secure storage
    of the work

-   **git checkout -b newbranchname –track origin/remotebranch**  
    Check out a branch stored on the remote repository, and allow
    pushing and pulling from the remote repository for that branch.

-   **git push**  
    Push the current branch to the remote branch that it is tracking.

-   **git merge master**  
    Merge the recent changes in master into your local branch (requires
    a pull first)

-   **git checkout master**  
    Switch branches to the master branch

-   **git merge newbranch**  
    Merge your code in ‘newbranch’ into your current branch (potentially
    master)

Each commit should contain one logical idea and the commit message
should clearly describe *everything* that is done in that commit. It is
fine for one commit to only contain a very minor change. Try and commit
regularly and avoid large commits. It is also a good idea to make sure
that code compiles before commiting. This helps catch errors that you
may be introducing and also allows the use of debugging tools such as
git bisect.

It should be noted that the ‘stable’ branch of the code, automatically
merged into from master upon successful completion of nightly tests, is
hosted on github on a public repository, and also pushed to the molpro
source code. The molpro developers will quickly send us angry emails if
poor code gets pushed into it from NECI, and I will be sure to forward
complaints onto the relevant parties!

