---
title: Tooling
---


### Using CTAGS with VIM

It is useful, especially for new developers, to be able to easily
navigate through NECI code. A simple solution for Vim users is to
generate a `tags` file containing the names of all functions and global
variables and their locations. Vim automatically reads this file from
the current directory, if it exists, and use it to facilitate code
navigation. Then you can jump to the deceleration of a variable or a
function by putting the cursor over it and pressing `Ctrl+]` . To go
back, press `Ctrl+t`. Other tag-related commands are explained here:
<https://vim.fandom.com/wiki/Browsing_programs_with_tags>.

To generate the tags file, a program called `ctags` is needed. It is
installed by default in many Linux distributions, but this version is
most probably the one called `ctags (GNU Emacs)` and does not support
the options we need. The required version is `Exuberant Ctags` or its
derivative `Universal Ctags` which you can be downloaded from here:
<https://github.com/universal-ctags/ctags>.

Once you installed the correct version and made sure it is the default
one, [^ctags] go to the source directory of NECI and run the script
`gen_vim_tags.sh` which is available in the tools directory

```bash
cd neci/src
../tools/gen_vim_tags.sh
```

This script does some tricks using the preprocessor to solve issues with
handling macros in NECI files. Without these, `ctags` would miss many
variables due to parsing issues. The script generates the necessary
`tags` file in the current directory and code navigation should become
available for all source files in this directory.

##### Note

`tags` is simply a text file with the symbol’s name, the file where its
defined, and the line number . There is no automatic magic happening
behind the scenes as one would expect from a full-fledged IDE.
*Therefore, whenever the code changes, you need to explicitly
re-generate the `tags` file*. Otherwise, Vim will simply jump to the old
positions of the symbols.

##### Tip

Another useful tool is `Tagbar` plugin for Vim which lists all
functions/variables in the current file in a side window. Using this
plugin does not require the `tags` file, because it generates its own
tags on-the-fly. However, you will need to install an updated version of
`ctags` executable. All necessary details are explained on the plugin’s
homepage: <https://github.com/majutsushi/tagbar>

##### EMACS Usage

Like Vim, Emacs accepts `TAGS` file (notice the capital letters). This
file can similarly be generated using `gen_emacs_tags.sh`, then you can
use `Meta+.` to jump to definition and `Meta+*` to go back. Other
tag-related commands are explained here:
<https://www.emacswiki.org/emacs/EmacsTags>.

