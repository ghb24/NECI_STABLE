#!/usr/bin/env python3

from collections import defaultdict
from os import walk
from os.path import join, abspath
import re
import sys


class TabCharacter:
    @staticmethod
    def test(line):
        return '\t' in line

    @staticmethod
    def report(file_path, line_number):
        return f'Tab characters in {file_path}, line {line_number}'

    @staticmethod
    def defined_for(file_path):
        return file_path.endswith(('.f', '.F', '.f90', '.F90', '.c'))


class TrailingSpace:
    blanks = re.compile(r'[ \t\r\f\v]+$')

    def test(self, line):
        return bool(self.blanks.search(line.strip('\n')))

    @staticmethod
    def report(file_path, line_number):
        return f'Trailing blanks in {file_path}, line {line_number}'

    @staticmethod
    def defined_for(file_path):
        return True


class WriteStar:
    write_star = re.compile(r'^ .*write\s*\(\s*\*', re.IGNORECASE)

    def test(self, line):
        return self.write_star.search(line)

    @staticmethod
    def report(file_path, line_number):
        return f'"write(*" in {file_path}, line {line_number}'

    @staticmethod
    def defined_for(file_path):
        return file_path.endswith(('.f', '.F', '.f90', '.F90'))


def run_tests_per_file(file_path, style_errors):
    errors = defaultdict(list)
    relevant = [x for x in style_errors if x.defined_for(file_path)]
    with open(file_path, 'r') as f:
        for line_number, line in enumerate(f):
            for style_error in relevant:
                if style_error.test(line):
                    errors[line_number].append(style_error)
    return dict(errors)


def run_tests(files, style_errors):
    errors = {}
    for file_path in files:
        errors[file_path] = run_tests_per_file(file_path, style_errors)
    return errors


def output_errors(errors):
    for file_path in errors:
        for line_number, style_errors_on_line in errors[file_path].items():
            for style_error in style_errors_on_line:
                print(style_error.report(file_path, line_number))


def get_files(start_dir):
    return (abspath(join(dname, fname)) for dname, dirs, fnames
            in walk(start_dir) for fname in fnames)


if __name__ == '__main__':
    STYLE_ERRORS = [TabCharacter(), TrailingSpace(), WriteStar()]
    files = get_files('../src')
    errors = run_tests(files, STYLE_ERRORS)
    output_errors(errors)
    if errors:
        sys.exit(1)
    else:
        sys.exit(0)

