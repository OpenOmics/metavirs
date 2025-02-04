#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Author: Skyler Kuhn

# Standard Library
from __future__ import print_function
import sys, os
import textwrap

# 3rd party packages,
# installed from pypi
import argparse
import pandas as pd


_help = textwrap.dedent("""create_matrix.py: Wake up Neo...
@Usage:
    $ ./create_matrix.py [-h] [--version] \\
            [--clean-suffix CLEAN_SUFFIX] \\
            [--nan-values NAN_VALUES] \\
            --join-on JOIN_ON \\
            --extract EXTRACT \\
            --output OUTPUT \\
            --input FILE_1 [FILE_2 ...]

@About:
    Given a list of files, an output file name, 
    a column name to join on, and a column name 
    containing features to extract, this script 
    will create a counts matrix for N sample by 
    M features. By default, sample column names
    in the resulting matrix will be the basename 
    of the X file that was provided; however, you 
    can provide string to strip any suffices from 
    the sample's filename.  

@Required Arguments:
    --input FILE_1 [FILE_2 ...]
                   One or more files to create the 
                   matrix from. Multiple files are 
                   joined on the `--join-on` column.
    --output OUTPUT  
                   Output file name. The file name of 
                   the resulting matrix file. Please
                   note that is resulting output file
                   will be tab-delimited. 
    --join-on JOIN_ON
                   Joins multiple files on this column
                   name. Please note this column name 
                   must exist in each `--input` file.  
    --extract EXTRACT 
                   Extracts features from this column 
                   to create the resulting matrix.

@Options
    --clean-suffix CLEAN_SUFFIX
                    Cleans each sample name in matrix to 
                    remove this suffix.
    --nan-values NAN_VALUES
                    Set missing feature values to this 
                    instead of the default. By default,
                    any missing fields are set to an
                    empty string, i.e ''. 
    --h, --help     Shows this help message and exits.


'The Matrix is everywhere. It is all around us. Even 
now in this very room' 
    - Morpheus 
"""
)


def err(*message, **kwargs):
    """Prints any provided args to standard error.
    kwargs can be provided to modify print functions 
    behavior.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    print(*message, file=sys.stderr, **kwargs)



def fatal(*message, **kwargs):
    """Prints any provided args to standard error
    and exits with an exit code of 1.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    err(*message, **kwargs)
    sys.exit(1)


def build_matrix(file, matrix_dict, join_index, parse_index, clean_suffix=None):
    """Creates a 2-D matrix, represented as a dictionary, where
    each key is N sample and its value is a nested dictionary 
    containing M feature and its extracted value.
    @param file <str>: 
        Path to input file with values to extract via `parse_index`.
    @param matrix_dict <dict[str][str]=[str]>: 
        Matrix dictionary to create where:
            [sample][feature] = feature_value
    @param join_index <int>: 
        Index of the field to join multiple files, usually the first column.
    @param parse_index <int>: 
        Index of field of interest to extract, could be counts.
    @param clean_suffix <str>:
        A suffix to remove from a sample name, could be a file extension.
        Default: None
	"""
    # Parse the input file to extract 
    # M features of interest for N sample
    with open(file, 'r') as fh:
        header = next(fh)               # skip over header 
        sample = os.path.basename(file) # remove path from sample
        if clean_suffix:
            # maybe replace with regex later
            sample = sample.split(clean_suffix)[0]   
        for line in fh:
            linelist = line.lstrip().rstrip().split('\t')
            j_field = linelist[join_index]
            j_value = linelist[parse_index]
            if sample not in matrix_dict:
                matrix_dict[sample] = {}
            matrix_dict[sample][j_field] = j_value
    return matrix_dict


def main():
    """Collect command line args and build the matrix."""
	# Parse command-line arguments
    # Create a top-level parser
    parser = argparse.ArgumentParser(
        usage = argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = _help,
        add_help=False
    )

    # Required Positional Arguments
    # List of input files
    parser.add_argument(
        '--input',
        required=True,
        nargs = '+',
        help = argparse.SUPPRESS
    )
    # Output file name
    parser.add_argument(
        '--output',
        required=True,
        type=str,
        help = argparse.SUPPRESS
    )
    # Join-on column name
    parser.add_argument(
        '--join-on',
        required=True,
        type=str,
        help = argparse.SUPPRESS
    )
    # Extract column name
    parser.add_argument(
        '--extract',
        required=True,
        type=str,
        help = argparse.SUPPRESS
    )

    # Options
    # Adding verison information
    parser.add_argument(
        '--clean-suffix', 
        required=False,
        default=None,
        type=str,
        help = argparse.SUPPRESS
    )

    # Set missing values to this
    # value instead of an empty 
    # string 
    parser.add_argument(
        '--nan-values', 
        required=False,
        default='',
        type=str,
        help = argparse.SUPPRESS
    )

    # Add custom help message
    parser.add_argument(
        '-h', '--help', 
        action='help', 
        help=argparse.SUPPRESS
    )

    # Collect parsed arguments
    args = parser.parse_args()

    # Sanity check for usage
    if len(sys.argv) == 1:
        # Nothing was provided
        fatal('Invalid usage: create_matrix.py [-h] ...')

    # Parse N sample's column 
    # of interest for M features  
    matrix = {}
    # Sort input file names to 
    # create more deterministic 
    # output regardless of the 
    # input order the files  
    for file in sorted(args.input):
        # Get first line of file
        # to get column indices 
        with open(file, 'r') as fh:
            header = fh.readline().lstrip().rstrip().split('\t')
            try:
                extract_field = header.index(args.extract)
                join_on = header.index(args.join_on)
            except ValueError:
                err(
                    'Error: Missing one of the follow fields: {} or {}'.format(
                        args.extract,
                        args.join_on
                    )
                )
                err('in the following --input {} file'.format(file))
                fatal('Please check that each file contains each column name!')

        # Create a matrix of N samples 
        # by M features 
        matrix = build_matrix(
            file=file, 
            matrix_dict=matrix,
            join_index=join_on,
            parse_index=extract_field,
            clean_suffix=args.clean_suffix
        )

    
    # Write matrix to output file
    df = pd.DataFrame(matrix)
    df.to_csv(
        args.output, 
        sep="\t", 
        header=True, 
        index=True, 
        index_label = args.join_on,
        na_rep=args.nan_values
    )


if __name__ == '__main__':
    # Create the matrix
    try:
        main()
    except BrokenPipeError:
        pass