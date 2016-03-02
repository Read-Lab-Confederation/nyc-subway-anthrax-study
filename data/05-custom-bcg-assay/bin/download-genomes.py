#! /usr/bin/env python

"""
Author: Robert A Petit III.

Date: 01/06/2016

Query NCBI Nucleotide database and retrieve results individually.
"""

import os
import time
import argparse as ap
from Bio import Entrez
import sys
Entrez.email = 'robert.petit@emory.edu'

if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='download-genomes.py', conflict_handler='resolve',
        description="Download sequences individually."
    )

    parser.add_argument('query', type=str, metavar="QUERY",
                        help=('Query to search.'))
    parser.add_argument('output', type=str, metavar="OUTPUT",
                        help=('File to write sequences to.'))
    parser.add_argument('--retmax', default=1000, type=int, metavar="INT",
                        help=('Maximum number of genomes to download '
                              '(Default: 1000)'))
    parser.add_argument('--skip', type=str,
                        help=('List of accessions to skip.'))
    parser.add_argument('--dry_run', action="store_true",
                        help=('Run as normal, but do not download any data.'))

    args = parser.parse_args()

    db = 'nuccore'
    retmax = 1000 if args.retmax == 0 else args.retmax
    term = args.query
    outdir = args.output
    rettype = 'fasta'
    retmode = 'text'

    # Also allow a list of accessions to be downloaded
    if os.path.isfile(args.query):
        with open(args.query, 'r') as fh:
            accessions = []
            for line in fh:
                accessions.append(line.rstrip())
        term = ('({0}) AND refseq[filter]'.format(' OR '.join(accessions)))

    skip_me = []
    if args.skip:
        with open(args.skip, 'r') as fh:
            for line in fh:
                skip_me.append(line.rstrip())

    handle = Entrez.esearch(db=db, retmax=retmax, term=term)
    esearch = Entrez.read(handle)
    if not args.dry_run:
        print "Database: {0}".format(db)
        print "Max Records: {0}".format(retmax)
        print "Output Directory: {0}".format(outdir)
        print "Query: {0}".format(term)
        print "----------"
        print "Searching for records..."
        print "\tFound {0} records.\n".format(esearch["Count"])
        print "Downloading records..."

    accessions = []
    for uuid in esearch['IdList']:
        handle = Entrez.esummary(db=db, id=uuid)
        esummary = Entrez.read(handle)

        if esummary[0]["Caption"] not in skip_me:
            accessions.append(esummary[0]["Caption"])
            if args.dry_run:
                # Do not download just print info
                print '{0}\t{1}'.format(
                    esummary[0]["Caption"],
                    esummary[0]["Title"].split(',')[0]
                )
            else:
                out_fa = '{0}/{1}.fasta'.format(outdir, esummary[0]["Caption"])
                if not os.path.isfile(out_fa):
                    print '\tDownloading {0}'.format(esummary[0]["Caption"])
                    efetch = Entrez.efetch(db=db, id=uuid, rettype=rettype,
                                           retmode=retmode)
                    with open(out_fa, 'w') as fh:
                        fh.write(efetch.read())
                    efetch.close()
                    # Don't get office banned! (again...)
                    print '\tSleeping for 5 seconds...'
                    time.sleep(5)
                else:
                    print '\tSkip existing {0}'.format(esummary[0]["Caption"])
        else:
            print "{0} is in the skip list".format(esummary[0]["Caption"])

    if not args.dry_run:
        print "Outputting list of completed genomes..."
        output = '{0}/completed-genomes.txt'.format(outdir)
        with open(output, 'w') as fh:
            fh.write('\n'.join(accessions))
