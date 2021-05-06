#!/bin/bash

# replace html colour commands with latex colour commands

# TODO I'm not sure if there is a better way to do this without just copying it
tmpfile=tmp_reformat_$1
cp $1 $tmpfile

sed -i 's/\[TOC\]/\\newpage/g' $tmpfile
# sed -i 's/@end\([^ ]*\)/\\End\{\1\}/g' $tmpfile
# sed -i 's/@\([^ ]*\)/\\Begin\{\1\}/g' $tmpfile # important that end is first
# spent way too long figuring these out only to realise it messes up the emails...
# guess I will hard code them.
# sed -i 's/@end\([^ ]*\)/\\dummy\{\\End\{\1\}\}/g' $tmpfile
# sed -i 's/@\([^ ]*\)/\\dummy\{\\Begin\{\1\}\}/g' $tmpfile # important that end is first
sed -i 's/@endnote/\\dummy\{\\End\{note\}\}/g' $tmpfile
sed -i 's/@note/\\dummy\{\\Begin\{note\}\}/g' $tmpfile 
sed -i 's/@endwarning/\\dummy\{\\End\{warning\}\}/g' $tmpfile
sed -i 's/@warning/\\dummy\{\\Begin\{warning\}\}/g' $tmpfile 
sed -i 's/@endtodo/\\dummy\{\\End\{todo\}\}/g' $tmpfile
sed -i 's/@todo/\\dummy\{\\Begin\{todo\}\}/g' $tmpfile 
sed -i 's/@endbug/\\dummy\{\\End\{bug\}\}/g' $tmpfile
sed -i 's/@bug/\\dummy\{\\Begin\{bug\}\}/g' $tmpfile 
sed -i 's/<span style="color: \([^<]*\)">\([^<]*\)<\/span>/\\textcolor{\1}{\2}/g' $tmpfile
# sed -i 's/(@[a-z]*)+/\\begin\{$1\}/g' $tmpfile

# CONVERT

# pandoc "$1" \
#        -f markdown-tex_math_dollars+tex_math_single_backslash \ 
#        -V geometry:margin=3cm \
#        -V geometry:a4paper \ 
#        -V fontsize=18pt \
#        -V numbersections=true \
#        --toc \
#        -o "$2" \
#        -t pdf 

# no idea why this works but the above doesn't (found online)
# pandoc "$1" \
#     -f markdown-tex_math_dollars+tex_math_single_backslash \
#     -V linkcolor:blue \
#     -V geometry:a4paper \
#     -V geometry:margin=2cm \
#     -V mainfont="DejaVu Serif" \
#     -V monofont="DejaVu Sans Mono" \
#     --pdf-engine=xelatex \
#     -o "$2"

# this works so clearly something is wrong in how I separate a multiline command
# TODO make this multiline...
pandoc -V geometry:margin=3cm -V geometry:a4paper -V fontsize=18pt -V numbersections=true --toc --pdf-engine=xelatex --include-in-header packages_to_include.tex "$tmpfile" -o "$2" -f markdown-tex_math_dollars+tex_math_single_backslash -t pdf 

# rm $tmpfile