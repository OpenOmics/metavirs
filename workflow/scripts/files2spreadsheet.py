#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Standard library
from __future__ import print_function
import argparse, sys, os, textwrap


# 3rd party packages,
# installed from pypi
import pandas as pd
import xlsxwriter

# Script metadata
__author__  = "Skyler Kuhn"
__version__ = "v0.2.0"

_help = textwrap.dedent(
"""./file2spreadsheet.py:
Creates an excel file from a list of input files.

@Usage:
    $ ./file2spreadsheet.py [-h] [--version] \\
            [--add-auto-filters] \\
            [--rm-suffix RM_SUFFIX] \\
            [--comment-symbol COMMENT_SYMBOL] \\
            --input FILE_1 [FILE_2 ...] \\
            --output OUTPUT

@About:
    Given a list of input files and an output file 
    name, this script will create an excel spread-
    sheet containing a tab/worksheet of each input
    file. 
    
    By default, any lines starting with the comment 
    character, '#', will not be skipped. Optionally,
    a suffix string can be provided to remove any 
    file extensions from each of the resulting tab
    or worksheet names. By default, the file exten-
    sion will only be removed from the when creating
    the tab or worksheet name.

@Required Arguments:
    -i, --input FILE_1 [FILE_2 ...]
                   One or more files to include in
                   the excel spreadsheet. Multiple 
                   files should be seperated by a
                   white space character.
    -o, --output OUTPUT  
                   Output excel file name. The file 
                   name of the resulting excel file. 
                   Within the excel file, there will 
                   be a tab for each provided input 
                   file. 
@Options
    -r, --rm-suffix CLEAN_SUFFIX
                    Removes the provided suffix the
                    tab/worksheet name in the excel 
                    file. By default, just the first
                    extension is removed.
    -c, --comment-symbol COMMENT_SYMBOL
                    Set the character that represents
                    comment lines. Lines starting with
                    this character will be skipped. By
                    default, all lines of the file will
                    be included and nothing is skipped.

    -a, --add-auto-filters
                    Adds auto-filters to all columns 
                    in a worksheet. Excel auto-filters 
                    apply drop-down filters/selectors
                    to the column headers in a sheet.
                    
    --h, --help     Shows this help message and exits.


@Author: {0}
@Version: {1}
""".format(__author__, __version__)
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


def parse_arguments():
    """Collect and parse command line arguments."""
    # Sanity check for usage
    if len(sys.argv) == 1:
        # Nothing was provided
        fatal('Invalid usage: file2spreadsheet.py [-h] ...')

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
        '-i',
        '--input',
        required=True,
        nargs = '+',
        help = argparse.SUPPRESS
    )
    # Output file name
    parser.add_argument(
        '-o',
        '--output',
        required=True,
        type=str,
        help = argparse.SUPPRESS
    )

    # Options
    # Skip lines with comment
    # chars, i.e "#" symbols
    parser.add_argument(
        '-c',
        '--comment-symbol', 
        required=False,
        default = None,
        type=str,
        help = argparse.SUPPRESS
    )

    # Remove suffix from 
    # tab/worksheet name
    parser.add_argument(
        '-r',
        '--rm-suffix', 
        required=False,
        default='',
        type=str,
        help = argparse.SUPPRESS
    )

    # Applies auto-filters  
    # to column headers
    parser.add_argument(
        '-a',
        '--add-auto-filters', 
        required=False,
        default=False,
        action = 'store_true',
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

    return args 


def reader(filename, subset=[], skip='#', **kwargs):
    """Reads in an MAF-like file as a dataframe. Determines the 
    correct handler for reading in a given MAF file. Supports reading
    in TSV files (.tsv, .txt, .text, .vcf, or .maf), CSV files (.csv), 
    and excel files (.xls, .xlsx, .xlsm, .xlsb, .odf, .ods, .odt ). 
    The subset option allows a users to only select a few columns 
    given a list of column names.
    @param filename <str>:
        Path of an MAF-like file to read and parse
    @param subset list[<str>]:
        List of column names which can be used to subset the df
    @param skip <str>:
        Skips over line starting with this character
    @params kwargs <read_excel()>
        Key words to modify pandas.read_excel() function behavior
    @return <pandas dataframe>:
        dataframe with spreadsheet contents
    """
    # Get file extension
    extension = os.path.splitext(filename)[-1].lower()

    # Assign a handler to read in the file
    if extension in ['.xls', '.xlsx', '.xlsm', '.xlsb', '.odf', '.ods', '.odt']:
        # Read in as an excel file
        return excel(filename, subset, skip, **kwargs)
    elif extension in ['.csv']:
        # Read in as an CSV file
        return csv(filename, subset, skip, **kwargs)
    else:
        # Default to reading in as an TSV file
        # Tab is the normal delimeter for MAF or VCF files
        # MAF files usually have one of the following
        # extensions: '.tsv', '.txt', '.text', '.vcf', '.maf'
        return tsv(filename, subset, skip, **kwargs)


def excel(filename, subset=[], skip='#', **kwargs):
    """Reads in an excel file as a dataframe. The subset option
    allows a users to only select a few columns given a list of 
    column names.
    @param filename <str>:
        Path of an EXCEL file to read and parse
    @param subset list[<str>]:
        List of column names which can be used to subset the df
    @param skip <str>:
        Skips over line starting with this character
    @params kwargs <read_excel()>
        Key words to modify pandas.read_excel() function behavior
    @return <pandas dataframe>:
        dataframe with spreadsheet contents
    """
    if subset:
        return pd.read_excel(filename, comment=skip, **kwargs)[subset]

    return pd.read_excel(filename, comment=skip, **kwargs)


def tsv(filename, subset=[], skip='#', **kwargs):
    """Reads in an TSV file as a dataframe. The subset option
    allows a users to only select a few columns given a list of 
    column names.
    @param filename <str>:
        Path of an TSV file to read and parse
    @param subset list[<str>]:
        List of column names which can be used to subset the df
    @param skip <str>:
        Skips over line starting with this character
    @params kwargs <read_excel()>
        Key words to modify pandas.read_excel() function behavior
    @return <pandas dataframe>:
        dataframe with spreadsheet contents
    """
    if subset:
        return pd.read_table(filename, comment=skip, **kwargs)[subset]

    return pd.read_table(filename, comment=skip, **kwargs)


def csv(filename, subset=[], skip='#', **kwargs):
    """Reads in an CSV file as a dataframe. The subset option
    allows a users to only select a few columns given a list of 
    column names.
    @param filename <str>:
         Path of an CSV file to read and parse
    @param subset list[<str>]:
        List of column names which can be used to subset the df
    @param skip <str>:
        Skips over line starting with this character
    @params kwargs <read_excel()>
        Key words to modify pandas.read_excel() function behavior
    @return <pandas dataframe>:
        dataframe with spreadsheet contents
    """
    if subset:
        return pd.read_csv(filename, comment=skip, **kwargs)[subset]

    return pd.read_csv(filename, comment=skip, **kwargs)


def excel_writer(files, spreadsheet,  skip_comments=None, remove_suffix = '', add_auto_filters=False):
    """Takes a list of files and creates one excel spreadsheet.
    Each file will becomes a sheet in the spreadsheet where the 
    name of the sheet is the basename of the file with the extension
    removed.
    @param files list[<str>]:
        List of files to merge into one execl file
    @param spreadsheet <str>:
        Output filename of the spreadsheet
    """
    # Create output directory as needed
    outdir = os.path.dirname(os.path.abspath(spreadsheet))
    if not os.path.exists(outdir):
        # Pipeline output directory
        # does not exist on filesystem
        os.makedirs(outdir)
    
    # Create a spreadsheet from the contents of each file
    with pd.ExcelWriter(spreadsheet, engine='xlsxwriter') as writer:
        for file in files:
            print('Reading in {}'.format(file))
            df = reader(file, skip=skip_comments)
            # Removing leading # char from header
            # df.rename(columns={df.columns[0]: "{0}".format(df.columns[0].lstrip('#'))})
            sheet, ext = os.path.splitext(os.path.basename(file))
            if remove_suffix:
                # Default behavior is to always
                # remove the file extension, how-
                # ever, if a user provides an 
                # additional suffix to remove
                # we must take into account the
                # extension has already been
                # removed, see above
                cleaned_suffix, ext = os.path.splitext(remove_suffix)
                sheet = sheet.rsplit(cleaned_suffix, 1)[0]
            
            # Sheet name cannot exceed 
            # 31 characters in length
            df.to_excel(writer, sheet_name = sheet[:31], index = False, freeze_panes = (1,0))
            # Automagically set the width
            # of a column to its max width
            worksheet = writer.sheets[sheet]  # pull worksheet object
            for idx, col in enumerate(df):  # loop through all columns
                series = df[col]
                max_len = max(
                    (
                        series.astype(str).map(len).max(),  # len of largest item
                        len(str(series.name))  # len of column name/header
                    )
                ) + 2  # adding a little extra space
                worksheet.set_column(idx, idx, max_len)  # set column width
                if add_auto_filters:
                    worksheet.autofilter(0, 0, df.shape[0], df.shape[1])


def main():
    # Parse command-line arguments
    args = parse_arguments()

    # List of files to convert into 
    # an excel spreasheet
    inputs = args.input
    # Output file name 
    output = args.output

    # Skip comment lines
    comment_char = args.comment_symbol
    # Remove suffix from 
    # X tab/worksheet name 
    rm_suffix = args.rm_suffix

    # Apply auto filters to 
    # column headers in a sheet
    auto_filter = args.add_auto_filters

    # Create XLSX file from the list
    # of input files 
    excel_writer(
        files=inputs,
        spreadsheet=output,
        skip_comments=comment_char,
        remove_suffix=rm_suffix,
        add_auto_filters=auto_filter
    )


if __name__ == '__main__':
    # Call main method
    main()
