#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Author: Skyler Kuhn

# Standard Library
from __future__ import print_function
import sys, os
import textwrap

# 3rd party packages,
# Install from pypi
import argparse # added in python/3.0


_help = textwrap.dedent("""collapse_on_taxa.py:
@Usage:
    $ ./collapse_on_taxa.py [-h] [--version] \\
            --collapse-on-taxa {family,genus,species} \\
            --output OUTPUT \\
            --input FILE

@About:
                        
    Given an annotated viral table file, an output
    file name, and a taxa column name to collapse 
    on, this script aggregate counts and coverage
    information for a given sample.

@Required Arguments:
    --input FILE
                   An annotated input viral table file.
                   This file should contain column the
                   following columns: 
                        • contig
                        • cov
                        • count
                        • taxid
                        • genus, family, or species
                   Taxon IDs were annotated with taxonkit.
    --output OUTPUT  
                   Output file name. The file name of 
                   the resulting matrix file. Please
                   note that is resulting output file
                   will be tab-delimited. 
    --collapse-on-taxa family,genus,species
                   Taxa column name to collapse info on.
                   If more than one contig is annotated 
                   to the same taxonomy, its counts and 
                   cov will be summated. Please note the
                   value provided to this option must 
                   exist as a column in the `--input` 
                   file.

@Options
    --h, --help     Shows this help message and exits.
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


def index_taxa(file, collapse_on_taxa, collapse_features = ['contig', 'count', 'cov', 'taxid']):
    """Indexes taxa information into a dictionary, where
    each key is a taxa name and its value is a nested dictionary 
    containing 1:M information of the following fields: 
    contig, count, cov, taxid
    @param file <str>: 
        Path to input file with values to extract/collapse via `collapse_on_taxa`.
    @param collapse_on_taxa <str>:
        Column name to collapse taxanomic information on.
        Examples: family or genus or species, etc
    @params collapse_features list[<str>]:
        Features to aggregate for collapsing.
        Default: ['contig', 'count', 'cov', 'taxid']
    @return collapsed_dict dict[`taxa_name`][`collapse_feature`] = [val1,val2,...]
        Example: 
        {'Betabaculovirus': 
            {'contig': ['k141_1852_flag=0_multi=1.0000_len=316', 'k142_1853_flag=0_multi=1.0000_len=316'], 
            'count': ['63', '123'], 
            'cov': ['1.0000', '2.0000'], 
            'taxid': ['56947', '56947']}, ...}
	"""
    # Parse the input file features
    # of interest
    collapsed_dict = {}
    with open(file, 'r') as fh:
        header = fh.readline().lstrip().rstrip().split('\t') 
        taxa_idx = header.index(collapse_on_taxa)
        for line in fh:
            linelist = line.lstrip().rstrip().split('\t')
            # Taxa name to collapse info on
            taxa_name = linelist[taxa_idx]
            if taxa_name not in collapsed_dict:
                # Initialize nested dict
                collapsed_dict[taxa_name] = {}
                for f in collapse_features:
                    collapsed_dict[taxa_name][f] = []
            # Append collapsed feature values to list
            for f in collapse_features:
                collapsed_dict[taxa_name][f].append(linelist[header.index(f)])
    return collapsed_dict


def main():
    """Collect command line args and collape info per-sample."""
	# Parse command-line arguments
    # Create a top-level parser
    parser = argparse.ArgumentParser(
        usage = argparse.SUPPRESS,
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = _help,
        add_help = False
    )

    # Required Positional Arguments
    # List of input files
    parser.add_argument(
        '--input',
        required = True,
        type = str,
        help = argparse.SUPPRESS
    )
    # Output file name
    parser.add_argument(
        '--output',
        required = True,
        type = str,
        help = argparse.SUPPRESS
    )
    # Collapse on taxa column name
    parser.add_argument(
        '--collapse-on-taxa',
        required = True,
        type = str,
        help = argparse.SUPPRESS,
        choices = [
            'kingdom',
            'phylum',
            'class',
            'order',
            'family',
            'genus',
            'species',
            'strain'
        ] 
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
        fatal('Invalid usage: collapse_on_taxa.py [-h] ...')

    # Get indices of relevent columns
    # for later parsing
    file = args.input
    with open(file, 'r') as fh:
        header = fh.readline().lstrip().rstrip().split('\t')
    
    # Sanity check to see if required 
    # columns are present in input file
    collapse_features = ["cov", "count", "taxid", "contig"]
    required_columns = [args.collapse_on_taxa] + collapse_features
    try:
        for c in required_columns:
            header.index(c)
    except ValueError:
        # Mising a required column
        err('Error: Missing one of the follow columns: {0}'.format(c))
        err('in the following --input {} file'.format(file))
        fatal('Please check that each file contains each column name!')

    # Index taxa information
    # key is name within a given
    # taxa and its value is a dict
    # with collapsed concat(contig), 
    # sum(count), sum(cov), and 
    # concat(taxid)
    taxa_idx = index_taxa(
        file = file,
        collapse_on_taxa = args.collapse_on_taxa,
        collapse_features = collapse_features
    )

    with open(args.output, 'w') as ofh:
        # Add header to output file
        ofh.write("{0}\tncontig\n".format('\t'.join([args.collapse_on_taxa] + collapse_features)))
        for taxa_name, aggr_features in taxa_idx.items():
            collapsed_contigs = ','.join(aggr_features['contig'])        # concat string
            num_contigs = len(aggr_features['contig'])                   # count concat contigs
            collapsed_taxids  = ','.join(aggr_features['taxid'])         # concat string
            sum_count = sum([ int(c) for c in aggr_features['count'] ])  # summate counts
            sum_cov   = sum([ float(c) for c in aggr_features['cov'] ])  # summate cov values
            # Write collapse taxa information
            # to the output file
            ofh.write(
                "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                    taxa_name,
                    sum_cov,
                    sum_count,
                    collapsed_taxids,
                    collapsed_contigs,
                    num_contigs
                )
            )


if __name__ == '__main__':
    # Collapse information on the given taxa level
    try:
        main()
    except BrokenPipeError:
        err("WARNING: Caught except related to pipefail! Did you use control-C?")