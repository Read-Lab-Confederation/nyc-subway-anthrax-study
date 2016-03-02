#! /usr/bin/env python
"""
    Author: Robert A Petit III
    Date: 01/06/2016

    Search NCBI Nucleotide database for completed Bacillus cereus Group
    genomes and download the FASTA file for each genome.
"""
import time
import argparse as ap
from Bio import Entrez
from urllib2 import HTTPError
Entrez.email = 'robert.petit@emory.edu'

if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='entrez-batch-download.py', conflict_handler='resolve',
        description="Download non B. cereus Group nucleotide in batch."
    )

    parser.add_argument('query', type=str, metavar="QUERY",
                        help=('Query to search.'))
    parser.add_argument('output', type=str, metavar="OUTPUT",
                        help=('File to write sequences to.'))
    parser.add_argument('--retmax', default=100000, type=int, metavar="INT",
                        help=('Maximum number of seqeunces to download '
                              '(Default: 100000)'))
    parser.add_argument('--batch_size', default=500, type=int, metavar="INT",
                        help=('Maximum number of seqeunces to download at a'
                              'time. (Default: 500)'))

    args = parser.parse_args()
    db = 'nuccore'
    term = args.query
    retmax = args.retmax
    batch_size = args.batch_size
    rettype = 'fasta'
    retmode = 'text'

    print "Database: {0}".format(db)
    print "Max Records: {0}".format(retmax)
    print "Output File: {0}".format(args.output)
    print "Query: {0}".format(term)
    print "----------"

    print "Searching for records..."
    handle = Entrez.esearch(db=db, retmax=retmax, term=term, usehistory="y")
    esearch = Entrez.read(handle)
    handle.close()

    print "\tFound {0} records.\n".format(esearch["Count"])
    print "Downloading records..."
    with open(args.output, 'w') as fh:
        for start in range(0, int(esearch["Count"]), batch_size):
            end = min(int(esearch["Count"]), start + batch_size)
            print "\tDownloading records {0} to {1}...".format(start + 1, end)
            attempt = 1
            error = True
            while attempt <= 1000 and error:
                error = True
                try:
                    efetch = Entrez.efetch(
                        db=db, rettype=rettype, retmode=retmode,
                        retstart=start, retmax=batch_size,
                        webenv=esearch["WebEnv"],
                        query_key=esearch["QueryKey"]
                    )
                    data = efetch.read()
                    if 'ERROR' in data:
                        print '\t\tDownload failed, sleep for 5s and retry...'
                        attempt += 1
                    else:
                        error = False
                        print '\t\tWriting to file and sleeping for 5s...'
                        fh.write(data)
                    efetch.close()
                    time.sleep(5)
                except HTTPError as err:
                    if 500 <= err.code <= 599:
                        print("Received error from server %s" % err)
                        print("Attempt %i of 3" % attempt)
                        attempt += 1
                        time.sleep(15)
                    else:
                        raise
