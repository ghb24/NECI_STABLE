#!/bin/bash
rm -f unsorted_tags
find . \( -iname \*.F -o -iname \*.F90 \)  -execdir cpp -I . -I .. -w -D "HElement_t=real" {} {}.tmp \; -exec ctags -a  -f "unsorted_tags" --sort=no --line-directives=yes --language-force=fortran {}.tmp \; -exec rm {}.tmp \;
find . -iname \*.template -exec ctags -a  -f "unsorted_tags" --sort=no  --language-force=fortran {} \; 
(head -n 7 unsorted_tags && tail -n +8 unsorted_tags | sort) > tags
rm unsorted_tags
