#!/bin/bash

# TODO I'm not sure if there is a better way to do this without just copying it
tmpfile=tmp_reformat_$1
cp $1 $tmpfile


# replace html colour commands with latex colour commands

sed -i 's/\[TOC\]/\\newpage/g' $tmpfile
sed -i 's/@endnote/\\dummy\{\\end\{note\}\}/g' $tmpfile
sed -i 's/@note/\\dummy\{\\begin\{note\}\}/g' $tmpfile
sed -i 's/@endwarning/\\dummy\{\\end\{warning\}\}/g' $tmpfile
sed -i 's/@warning/\\dummy\{\\begin\{warning\}\}/g' $tmpfile
sed -i 's/@endtodo/\\dummy\{\\end\{todo\}\}/g' $tmpfile
sed -i 's/@todo/\\dummy\{\\begin\{todo\}\}/g' $tmpfile
sed -i 's/@endbug/\\dummy\{\\end\{bug\}\}/g' $tmpfile
sed -i 's/@bug/\\dummy\{\\begin\{bug\}\}/g' $tmpfile
sed -i 's/<span style="color: \([^<]*\)">\([^<]*\)<\/span>/\\textcolor{\1}{\2}/g' $tmpfile

# CONVERT

pandoc \
      -V geometry:margin=3cm \
      -V geometry:a4paper \
      -V fontsize=18pt \
      -V numbersections=true \
      --toc \
      --include-in-header packages_to_include.tex \
      --pdf-engine=xelatex \
      "$tmpfile" \
      -o "$2" \
      -f markdown-tex_math_dollars+tex_math_single_backslash

rm $tmpfile
