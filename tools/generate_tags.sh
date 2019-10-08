#!/bin/bash

# Remove the temporary file 'unsorted_tags' in case it is sill here.
rm -f unsorted_tags

# 'ctags' has trouble handling preprocessor directives in fortran files. 
# Therefore, we call the preprocessor 'cpp' explicitly for each fortran file,
# store the result in a temporary file and use ctags to append its tags
# to 'unsorted_tags' file.
find . \( -iname \*.F -o -iname \*.F90 \) \
    -execdir cpp -I . -I .. -w -D "HElement_t=real" {} {}.tmp \; \
    -exec ctags --excmd=number -a  -f "unsorted_tags" \
    --sort=no --line-directives=yes --language-force=fortran {}.tmp \; \
    -exec rm {}.tmp \;

# Suprisingly, 'ctags' works well with template files
find . -iname \*.template \
    -exec  ctags -a  -f "unsorted_tags" \
    --excmd=number  --sort=no --language-force=fortran {} \; 

# For some reason, using 'ctags' above with '--sort=yes' option leads to 
# problems in the output. Therefore, we generated 'unsorted_tags' and sort it 
# here ouselves.
(head -n 7 unsorted_tags && tail -n +8 unsorted_tags | sort) > tags

# Remove the temporary unsorted file
rm unsorted_tags
