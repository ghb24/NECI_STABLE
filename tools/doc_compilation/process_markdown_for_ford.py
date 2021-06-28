#!/usr/bin/env python3


import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.customization import convert_to_unicode
from pathlib import Path
from functools import cmp_to_key
import argparse
import re
import sys


def ref_string(entry):
    """
    takes bibtexparser entry and pretty prints it
    """
    to_print = ''
    if 'title' in entry:
        to_print += '*'+entry['title'].strip().replace('\n', ' ')+'.* '
    if 'author' in entry:
        to_print += entry['author'].strip().replace(' and ', ', ')+'; '
    if 'journal' or 'journaltitle' in entry:
        if 'journal' in entry:
            to_print += entry['journal'].strip() + ' -- '
        elif 'journaltitle' in entry:
            to_print += entry['journaltitle'].strip() + ' -- '
        if 'year' in entry:
            to_print += '(' + entry['year'] + ') '
        elif 'date' in entry:
            to_print += '(' + entry['date'].split('-')[0] + ') '
        if 'volume' in entry:
            to_print += entry['volume']+' '
        if 'pages' in entry:
            to_print += entry['pages']

    return to_print.strip()

# TODO the two functions below share a lot in common, probably good to modularise


def replace_keys_get_ref(line, bib_database, to_end='', n=1, reffile=''):
    """
    given a string and a bibtexparser database object, gives back a parsed
    string based on the database. Also returns the modified database as well as
    the references section (along with n which signifies order)
    """
    cite_regex = re.compile('\[\@(.*?)\]')
    keys_in_line = cite_regex.findall(line)
    line_parsed = line
    # to_end = ''
    for k in keys_in_line:
        if k in bib_database.entries_dict:
            if 'n' not in bib_database.entries_dict[k]:
                bib_database.entries_dict[k].update({'n': n})
                to_end += '\n\n<a id=\"%s\"></a>[%i] ' % (
                    k, n) + ref_string(bib_database.entries_dict[k])
                n += 1
            line_parsed = line_parsed.replace(
                '[@'+k+']', '<a href=\"%s#%s\">[%i]</a>' % (reffile, k, bib_database.entries_dict[k]['n']))

    return line_parsed, bib_database, to_end, n


def replace_footnotes(line, n=1, to_end=''):
    foot_regex = re.compile('footnote\{(.*?)\}')
    footnotes_all = foot_regex.findall(line)
    line_parsed = line
    # n=1
    for footnote in footnotes_all:
        to_end += '\n\n<a id=\"%i\"></a><sup>%i</sup>' % (n, n) + footnote
        line_parsed = line_parsed.replace(
            '\\footnote{'+footnote+'}', '<a href=\"#%i\"><sup>%i</sup></a>' % (n, n))
        n += 1
    return line_parsed, to_end, n


def preprocess_markdown_file(f, bib_database, reffile='', n=1, refs=''):
    # TODO: processed_file could become a generator returning lines.
    processed_file = []
    break_yaml = '---'
    after_preamble = False

    nfoot = 1
    to_end = ''
    line = f.readline()
    while line:
        if after_preamble:
            # processed_file.append(line)
            line_parsed, bib_database, refs, n = replace_keys_get_ref(
                line, bib_database, to_end=refs, n=n, reffile=reffile)
            line_parsed, to_end, nfoot = replace_footnotes(
                line_parsed, n=nfoot, to_end=to_end)
        else:
            if line.strip() == '---':
                after_preamble = True
            line_parsed = line
        processed_file.append(line_parsed)
        line = f.readline()

    # add footnotes to end if applicable
    if to_end:
        processed_file.append('\n\n***'+to_end)

    return processed_file, bib_database, refs, n


def process_dir(bibfile, md_files, out_dir, build_dir):
    outlitfile = out_dir / Path('{}.md'.format(bibfile.stem))
    outlitfilehtml = build_dir / Path('page') / Path(out_dir.name) / Path('{}.html'.format(bibfile.stem))

    out_dir.mkdir(parents=True, exist_ok=True)

    with open(bibfile, 'r') as bibtex_file:
        parser = BibTexParser()
        parser.customization = convert_to_unicode
        bib_data = bibtexparser.load(bibtex_file, parser=parser)


    n = 1
    refs = ''
    for fname in md_files:
        if fname.resolve() == (out_dir / fname.name).resolve():
            raise ValueError("Script would overwrite the input. Choose different out_dir.")

        with open(fname, 'r') as f, open(out_dir / fname.name, 'w') as fp:
            processed_lines, bib_data, refs, n = preprocess_markdown_file(
                f, bib_data, reffile=outlitfilehtml, n=n, refs=refs)

            for l in processed_lines:
                fp.write(l)

    with open(outlitfile, 'w') as outfile:
        outfile.write('title: References\n')
        outfile.write('---')
        outfile.write(refs)


def driver(bibfile, md_dirs, out_dir, build_dir):
    @cmp_to_key
    def compare(p1, p2):
        "Compare two paths lexicographically, except if the file is named `index.md`"
        if p1 == p2:
            return 0
        elif p1.name == 'index.md':
            return -1
        elif p2.name == 'index.md':
            return 1
        else:
            return -1 if p1 < p2 else 1

    for dir in md_dirs:
        process_dir(bibfile, sorted(dir.glob('*.md'), key=compare), out_dir / Path(dir.name), build_dir)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bibfile', type=Path, help='The BibTex File', required=True)
    parser.add_argument('--out_dir', type=Path, help='The Output directory', required=True)
    parser.add_argument('--md_dirs', type=Path, help='The directories with md files', nargs='+', required=True)
    parser.add_argument('--build_dir', type=Path, help='The directory where html files of the final documentation reside.', required=True)

    args = parser.parse_args()

    driver(args.bibfile, args.md_dirs, args.out_dir, args.build_dir)
