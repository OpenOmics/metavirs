#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
from concurrent.futures import ThreadPoolExecutor
import time, argparse, os, sys, shutil

# Third-party pypi packages
# pip install as needed
import requests

# Local imports
from utils import (
    fatal,
    exists,
    err,
    md5sum
)


# Constasts
__version__ = 'v0.1.0'
__author__  = 'Skyler Kuhn'

# Functions
def retry(times=5, exceptions=(Exception)):
    """
    Decorator to retry running a function. Retries the wrapped function 
    N times with an exponential backoff stragety. A tuple of Exceptions
    can be passed that trigger a retry attempt. When times is equal to 
    4 the back-off strategy will be {4, 16, 64, 256} seconds. Calls fatal
    if the function cannot be run within the defined number of times.
    @param times <int>: 
        The number of times to repeat the wrapped function,
        default: 5
    @param exceptions tuple(<Exception>):
        Tuple of Python Exceptions that will trigger a retry attempt,
        default: (Exception)
    @return <object func>:
        Calls fatal when func cannot be run after N times.
    """
    def decorator(func):
        def do(*args, **kwargs):
            # Begin the attempt loop and return if successful 
            attempt = 0
            delay = 1
            while attempt < times:
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    # Default expection, Exception
                    err('Function failed: {0}'.format(func.__name__))
                    err('\t@args: {0}'.format(args))
                    err('\t@kwargs: {0}'.format(kwargs))
                    err('\t@reason: {0}'.format(e))
                    # Increase backoff: 4, 16, 64, 256, 1024, 4096...
                    attempt += 1
                    delay = 4**attempt
                    err("Attempt: {0}/{1}".format(attempt, times))
                    err("Trying again in {0} seconds!\n".format(delay))
                    time.sleep(delay)
            # Could not succesfully run the given function 
            # in its alloted number of tries, exit with error 
            fatal('Fatal: failed after max retry attempts!')
        return do
    return decorator


@retry()
def download(session, url, output_file=None, md5_checksum=None):
    """
    Download a file locally with a requests Session object.
    @param session <object requests.Session>: 
        Session object for persistant cookies and connection pooling
    @param url <str>:
        https url of the file to download from server
    @param output_file <str>:
        Output file name of file to download locally
    @return tuple(output_file <str>, hasher.hexdigest <str>):
        Output file name of downloaded file
        MD5 checksum of the dowloaded file
    """
    local_filename = output_file
    if not output_file:
        local_filename = url.split('/')[-1]

    with session.get(url, stream=True) as response:
        # Raise an  HTTPError object if an error 
        # has occurred during the process
        response.raise_for_status()
        print("Downloading {0} from {1}".format(local_filename, url))
        with open(local_filename, 'wb') as ofh:
            shutil.copyfileobj(response.raw, ofh)
        print("Completed download of {0}".format(local_filename))

    # Check downloaded MD5 against 
    # expected MD5 checksum value  
    local_checksum = ''
    if md5_checksum:
        local_checksum = md5sum(local_filename)
        if local_checksum != md5_checksum:
            raise Exception(
                'Retrying {0} download... expected {1} MD5, recieved {2}'.format(
                    local_filename,
                    md5_checksum,
                    local_checksum
                )
            )

    return local_filename, local_checksum


def cli():
    """
    Parses command-line arguments when script is directly called.
    @return args <object argparse ArgumentParser.parse_args()>:
    """
    # Command-line argument parser
    parser = argparse.ArgumentParser(
        description = 'download.py: Download large resource bundles in parallel'
    )

    # Verison information
    parser.add_argument(
        '--version', 
        action='version', 
        version='%(prog)s {}'.format(__version__)
    )

    # Required Arguments
    # Input files to download over HTTPS,
    # resources bundle is in datashare
    parser.add_argument(
        '--input',
        # Check if the file exists and if it is readable
        required = True,
        nargs = '+',
        help = 'Required: Input file URLs to download over HTTPS.'
    )
    # Output path to store downloaded files,
    # will be created if does not exist
    parser.add_argument(
        '--output',
        type = lambda option: os.path.abspath(os.path.expanduser(option)),
        required = True,
        help = 'Required: Output directory to store downloaded files.'
    )

    # Options
    # MD5 checksums of each input file,
    # the number of checksums must match
    # the number of provided inputs URLs
    parser.add_argument(
        '--md5',
        required = False,
        default = [],
        nargs = '+',
        help = 'Option: MD5 checksums of input files.'
    )

    # Threads for concurrent downloads

    parser.add_argument(
        '--threads',
        type = int,
        default = 2,
        required = False,
        help = 'Option: Number of threads for concurrent downloads.'
    )

    # Dry-run to see what files will
    # be downloaded locally
    parser.add_argument(
        '--dry-run',
        action = 'store_true',
        required = False,
        default = False,
        help = 'Option: Does not execute anything, \
            displays what files will be downloaded.'
    )

    # Forces re-downloading of all files, 
    # by default files that DNE are downloaded.
    parser.add_argument(
        '--force',
        action = 'store_true',
        required = False,
        default = False,
        help = 'Option: Force re-download of files \
            that already exist locally.'
    )

    args = parser.parse_args()
    
    return args


def main(options):
    # Check that number of MD5s
    # that were provided is equal
    # to the number of input URLs
    input2md5 = {}
    if options.md5:
        if len(options.md5) != len(options.input):
            fatal(
                'Usage: error: Number of --input URLs '
                'does not equal the number of provided '
                '--md5 CHECKSUMS!'
            )
        # Dictionary to map input 
        # to expected MD5 checksum 
        input2md5 = { 
            options.input[i]:options.md5[i] \
                for i in range(len(options.input))
        }

    # Create output directory,
    # if it does not exist
    if not exists(options.output):
        os.makedirs(options.output)

    # Determine what files need
    # to be downloaded, if they
    # already exist then do not
    # re-download the file
    urls = []
    for url in options.input:
        local_file = os.path.join(options.output, url.split('/')[-1])
        if local_file and not exists(local_file):
            urls.append(url)

    # Force re-downloading of all
    # files, even if they already
    # exist locally on the FS
    if options.force:
        urls = options.input

    # Create list of output files,
    # but first remove any deplicates
    dedup = [] 
    for url in urls:
        if url not in dedup:
            dedup.append(url)
    urls = dedup
    output_files = []
    for url in urls:
        local_file = os.path.join(
            options.output, 
            url.split('/')[-1]
        )
        output_files.append(local_file)

    # Display what files would be 
    # downloaded and exit
    if options.dry_run:
        print('Dry-running download.py:\n')
        print('Listing files to download...')
        for f in output_files:
            print("  •", f)
        return

    # Number of threads to use
    # for concurrent downloads
    nt = int(options.threads)
    nsessions = len(urls)
    md5checksums = []
    for url in urls:
        if options.md5:
            md5checksums.append(input2md5[url])
    # Create a thread pool to execute calls 
    # asynchronously using context-switching,
    # this circumnavigates python's GIL. 
    # NOTE: Deadlocks can occur when the 
    # callable associated with a Future waits 
    # on the results of another Future. 
    with ThreadPoolExecutor(max_workers=nt) as nthreads:
        # Create a Session object, this allows 
        # for persists cookies that can be used 
        # across all requests made from Session 
        # instance. It will also use urllib3’s 
        # connection pooling, so if you’re making 
        # several requests to the same host, the 
        # underlying TCP connection will be 
        # reused, which can result in large 
        # performance gains. 
        with requests.Session() as session:
            results = nthreads.map(
                download, 
                [session]*nsessions, 
                urls,
                output_files,
                md5checksums
            )
            nthreads.shutdown(wait=True)
        
    # Gather results
    for result in results:
        print(result)


if __name__ == '__main__':
    args = cli()
    main(options = args)
