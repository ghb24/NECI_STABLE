#!/bin/bash

# TODO I'm not sure if there is a better way to do this without just copying it
tmpfile=tmp_reformat_$1
cp $1 $tmpfile


# replace html colour commands with latex colour commands

sed -i 's/\[TOC\]/\\newpage/g' $tmpfile
sed -i 's/@endnote/\\dummy\{\\End\{note\}\}/g' $tmpfile
sed -i 's/@note/\\dummy\{\\Begin\{note\}\}/g' $tmpfile
sed -i 's/@endwarning/\\dummy\{\\End\{warning\}\}/g' $tmpfile
sed -i 's/@warning/\\dummy\{\\Begin\{warning\}\}/g' $tmpfile
sed -i 's/@endtodo/\\dummy\{\\End\{todo\}\}/g' $tmpfile
sed -i 's/@todo/\\dummy\{\\Begin\{todo\}\}/g' $tmpfile
sed -i 's/@endbug/\\dummy\{\\End\{bug\}\}/g' $tmpfile
sed -i 's/@bug/\\dummy\{\\Begin\{bug\}\}/g' $tmpfile
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
