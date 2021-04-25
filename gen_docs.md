project: NECI
output_dir: ./doc
project_github: https://github.com/ghb24/NECI_STABLE/
summary: for FCIQMC
author: Alavi Group
github: 
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
graph: false
search: false 
include: ./src
         ./src/lib
         ./utils
         ./tools
parallel: 4
dbg: true
warn: true
page_dir: pages

#Summary 

Some preliminary documentation generated using Ford. Over here we would probably 
want to include some kind of summary of what NECI is. Note there is decent 
documentation on [how to use Ford on its Github](https://github.com/Fortran-FOSS-Programmers/ford/wiki).

Note we can also have links to Github, Gitlab, a logo above. 

@note 
At present, I have found 2 issues with Ford: one is that it [does not support 
multi-line strings with exclamation marks](https://github.com/Fortran-FOSS-Programmers/ford/issues/320) (this is just a bug), so I had to do 
a workaround for this. Another issue is that it crashes for this project, though 
I am not confident why. I have managed to fix it by some troubleshooting. In 
order to build with Ford, first install as instructed then replace the file 
`ford/templates/macros.html` with the `macros.html` founds in `docs`. 
For reference, the few lines added are marked with the comment
`<!-- NOTE this wasn't originally here (PHILIP)-->`
@endnote

In order to generate the docs, simply use `ford gen_docs.md`. The search index 
is currently the bottleneck for document generation. If you want to be able to 
search the docs, simply replace `search: false` with `search: true` above (this
will take a bit longer to generate, though).

@note 
As warnings are enabled, Ford will output many warnings at the moment. This is because we do not yet have all the code formatted in a way
that is understandable to the program (Ford is warning about undocumented variables, of which there are many). The default format is `!!` but it can be changed. It is worth going through [this page](https://github.com/Fortran-FOSS-Programmers/ford/wiki/Writing-Documentation).
@endnote

Ford has some useful notes you may recognise from using Doxygen, like 

@note
notes 
@endnote

@warning 
warning tags 
@endwarning

@todo
todo tags
@endtodo

@bug 
bug tags 
@endbug