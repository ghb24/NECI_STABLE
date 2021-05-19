#!/bin/bash

# Invoke with two arguments
# First argument: Markdown input file
# Second argument Latex Header
# Third argument: PDF output file

md_input=$1
latex_header=$2
pdf_output=$3

# TODO I'm not sure if there is a better way to do this without just copying it
tmpfile=/tmp/tmp_reformat_neci_doc.md
cp ${md_input} $tmpfile


# replaces markdown [TOC] command with a latex page break command
# (pandoc already creates a TOC and we prefer the TOC on its own page(s))
sed -i 's/\[TOC\]/\\newpage/g' $tmpfile

# replace @note ... @endnote with \begin{note}  \end{note}
# The trick is that contents of LaTeX commands are not parsed by pandoc.
# The verbatimLaTeX macro is just mirroring its output.
sed -i 's/@endnote/\\verbatimLaTeX\{\\end\{note\}\}/g' $tmpfile
sed -i 's/@note/\\verbatimLaTeX\{\\begin\{note\}\}/g' $tmpfile
sed -i 's/@endwarning/\\verbatimLaTeX\{\\end\{warning\}\}/g' $tmpfile
sed -i 's/@warning/\\verbatimLaTeX\{\\begin\{warning\}\}/g' $tmpfile
sed -i 's/@endtodo/\\verbatimLaTeX\{\\end\{todo\}\}/g' $tmpfile
sed -i 's/@todo/\\verbatimLaTeX\{\\begin\{todo\}\}/g' $tmpfile
sed -i 's/@endbug/\\verbatimLaTeX\{\\end\{bug\}\}/g' $tmpfile
sed -i 's/@bug/\\verbatimLaTeX\{\\begin\{bug\}\}/g' $tmpfile

# replace with literature.md file contents (ford command)
# couldn't figure out the general command, but hard-coding literature.md is fine
sed -i "s/{!docs\/literature.md!}/$(awk 1 ORS='\\n' literature.md)/g" $tmpfile 
# replace ford markdown newline with pandoc newline
sed -i 's/<br>$/\\/g' $tmpfile

# replace html colour commands with latex colour commands
sed -i 's/<span style="color: \([^<]*\)">\([^<]*\)<\/span>/\\textcolor{\1}{\2}/g' $tmpfile

# CONVERT
pandoc \
      -V geometry:margin=3cm \
      -V geometry:a4paper \
      -V fontsize=18pt \
      -V numbersections=true \
      --toc \
      --include-in-header ${latex_header} \
      --listings \
      "$tmpfile" \
      -o "${pdf_output}" \
      -f markdown-tex_math_dollars+tex_math_single_backslash

rm $tmpfile
