#!/usr/bin/env python

from __future__ import print_function

from collections import defaultdict
from os import walk
from os.path import join
import re
import sys
import argparse


class TabCharacter:
    @staticmethod
    def test(line):
        return '\t' in line

    @staticmethod
    def report(file_path, line_number):
        return 'Tab characters in {}:{}'.format(file_path, line_number)

    @staticmethod
    def defined_for(path):
        return _fortran_suffix(path) or _NECI_template_suffix(path)


class TrailingSpace:
    blanks = re.compile(r'[ \t\r\f\v]+$')

    def test(self, line):
        return self.blanks.search(line.strip('\n'))

    @staticmethod
    def report(file_path, line_number):
        return 'Trailing blanks in {}:{}'.format(file_path, line_number)

    @staticmethod
    def defined_for(path):
        return (_fortran_suffix(path) or _C_suffix(path)
                or _NECI_template_suffix(path))


class WriteStar:
    write_star = re.compile(r'^ .*write\s*\(\s*\*', re.IGNORECASE)

    def test(self, line):
        return self.write_star.search(line)

    @staticmethod
    def report(file_path, line_number):
        return '"write(*" in {}:{}'.format(file_path, line_number)

    @staticmethod
    def defined_for(path):
        return _fortran_suffix(path) or _NECI_template_suffix(path)


class CompileFlag:
    compile_flag = re.compile(r'(\bifdef\b|\bifndef).*?__.*?', re.IGNORECASE)

    def test(self, line):
        return self.compile_flag.search(line)

    @staticmethod
    def report(file_path, line_number):
        return '"Compile Flag with Double Underscore" in {}:{}'.format(file_path, line_number)

    @staticmethod
    def defined_for(path):
        return _fortran_suffix(path) or _NECI_template_suffix(path) or _C_suffix(path)


def run_tests_per_file(file_path, style_errors):
    errors = defaultdict(list)
    relevant = [x for x in style_errors if x.defined_for(file_path)]
    try:
        with open(file_path, 'r') as f:
            for line_number, line in enumerate(f):
                for style_error in relevant:
                    if style_error.test(line):
                        errors[line_number].append(style_error)
    except UnicodeDecodeError:
        pass
    return dict(errors)


def run_tests(files, style_errors):
    errors = {}
    for file_path in files:
        errors_per_file = run_tests_per_file(file_path, style_errors)
        if errors_per_file:
            errors[file_path] = errors_per_file
    return errors


def output_errors(errors):
    for file_path in errors:
        for line_number, style_errors_on_line in errors[file_path].items():
            for style_error in style_errors_on_line:
                # python is 0 indexed. Go to 1 indexed line number for output
                print(style_error.report(file_path, line_number + 1))


def get_files(start_dir):
    return (join(dname, fname) for dname, dirs, fnames
            in walk(start_dir) for fname in fnames)


def parse_args():
    parser = argparse.ArgumentParser(description='src dir')
    parser.add_argument('src_dir', type=str, help='src directory to check')
    args = parser.parse_args()
    return args.src_dir


def _fortran_suffix(path):
    return path.endswith(('.f', '.F', '.f90', '.F90', '.fpp'))


def _C_suffix(path):
    return path.endswith('.c')


def _NECI_template_suffix(path):
    return path.endswith(
        tuple(('{}.template'.format(suff)
              for suff in ('.f', '.F', '.f90', '.F90'))))


if __name__ == '__main__':
    STYLE_ERRORS = [TabCharacter(), WriteStar(), CompileFlag()]
    src_dir = parse_args()
    files = get_files(src_dir)
    errors = run_tests(files, STYLE_ERRORS)
    output_errors(errors)
    if errors:
        sys.exit(1)
    else:
        sys.exit(0)
