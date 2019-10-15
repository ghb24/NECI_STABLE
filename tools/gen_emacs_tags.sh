#!/bin/bash                                                                     
                                                                                
# Remove 'TAGS' file in case it is sill here.                                   
rm -f TAGS                                                                      

# 'ctags' has trouble handling preprocessor directives in fortran files. 
# Therefore, we call the preprocessor 'cpp' explicitly for each fortran file,
# store the result in a temporary file and use ctags to append its tags
# to 'TAGS' file.
find . \( -iname \*.F -o -iname \*.F90 \) \
    -execdir cpp -I . -I .. -w -D "HElement_t=COMPLEX" {} {}.tmp \; \
    -exec ctags --excmd=number -e -a  -f "TAGS" \
    --sort=no --line-directives=yes --language-force=fortran {}.tmp \; \
    -exec rm {}.tmp \;

# Suprisingly, 'ctags' works well with template files
find . -iname \*.template \
    -exec  ctags -e -a  -f "TAGS" \
    --excmd=number  --sort=no --language-force=fortran {} \; 

# During preprocesing, we replaced HElement_t macros with COMPLEX to avoid 
# problems in parsing the source files. Here we replace them back in the TAGS
# file.
sed -i 's/COMPLEX/HElement_t/g' TAGS
