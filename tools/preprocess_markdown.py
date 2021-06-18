# %%
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.customization import convert_to_unicode
import re
import sys

# %%


def ref_string(entry):
    """
    takes bibtexparser entry and pretty prints it
    """
    to_print = ''
    if 'title' in entry:
        to_print += '*'+entry['title'].strip().replace('\n', ' ')+'* '
    if 'author' in entry:
        to_print += entry['author'].strip().replace(' and ', ', ')+' '
    if 'journal' in entry:
        to_print += entry['journal'].strip() + ' -- '
        if 'year' in entry:
            to_print += '(' + entry['year'] + ') '
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
    processed_file = []
    break_yaml = '---'
    after_preamble = False

    # n = 1
    # refs = ''
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


# %%

if __name__ == "__main__":
    bibfile = 'lit.bib'
    outlitfile = 'literature.md'
    outlitfilehtml = '.'.join(outlitfile.rsplit('.')[:-1])+'.html'
    md_files = sys.argv[1:]

    with open(bibfile) as bibtex_file:
        parser = BibTexParser()
        parser.customization = convert_to_unicode
        bib_data = bibtexparser.load(bibtex_file, parser=parser)

    # print(md_files)
    n = 1
    refs = ''

    for fname in md_files:
        with open(fname, 'r') as f:
            processed_lines, bib_data, refs, n = preprocess_markdown_file(
                f, bib_data, reffile=outlitfilehtml, n=n, refs=refs)
        with open(fname + "_parsed", 'w') as fp:  # TODO this is obviously a bad solution
            for l in processed_lines:
                fp.write(l)

    # NOTE I am overwriting the file if it's already there!
    with open(outlitfile, 'w') as outfile:
        outfile.write('title: References\n')
        outfile.write('---')
        outfile.write(refs)
