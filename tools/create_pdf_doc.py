#!/usr/bin/env python3

import argparse
import sys
import re
from pathlib import Path


ENVIRONMENTS = {'note', 'warning', 'todo', 'bug'}
START_ENVIRONMENTS = {
    '@{}'.format(env): '\\verbatimLaTeX{{\\begin{{{}}}}}'.format(env) for env in ENVIRONMENTS}
END_ENVIRONMENTS = {
    '@end{}'.format(env): '\\verbatimLaTeX{{\\end{{{}}}}}'.format(env) for env in ENVIRONMENTS}


def process_markdown_file(f):
    processed_file = []
    could_be_title = False

    color_regex = re.compile(r'<span style="color: (.*)">(.*)<\/span>')

    line = f.readline()
    n = 0
    while line:
        line = line.replace(r'<br>', '\\')

        if n == 4:
            processed_file.append('\\newpage\n')
            processed_file.append('# Using NECI\n')

        stripped = line.strip()
        if stripped in START_ENVIRONMENTS.keys():
            processed_file.append(START_ENVIRONMENTS[stripped] + '\n')
        elif stripped in END_ENVIRONMENTS.keys():
            processed_file.append(END_ENVIRONMENTS[stripped] + '\n')
        elif color_regex.search(line):
            processed_file.append(color_regex.sub(
                r'\\textcolor{\1}{\2}', line))
        elif line == '---\n' and len(processed_file) >= 3:
            if 'title:' in processed_file[-1] and '---\n' == processed_file[-2]:
                processed_file.pop()
                processed_file.pop()
            else:
                processed_file.append(line)
        else:
            processed_file.append(line)

        line = f.readline()
        n += 1
    return processed_file


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin)

    args = parser.parse_args()

    if not args.input_file:
        sys.exit("Please provide an input file, or pipe it via stdin")

    processed_file = process_markdown_file(args.input_file)

    for line in processed_file:
        print(line, end='')


if __name__ == '__main__':
    main()
