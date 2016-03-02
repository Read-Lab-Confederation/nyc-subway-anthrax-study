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
    parser.add_argument('--skip', type=str, metavar="OUTPUT",
                        help=('List of accessions to skip.'))

    args = parser.parse_args()
    db = 'nuccore'
    retmax = 1000 if args.retmax == 0 else args.retmax
    term = args.query
    outdir = args.output
    rettype = 'fasta'
    retmode = 'text'

    skip_me = []
    if args.skip:
        with open(args.skip, 'r') as fh:
            for line in fh:
                skip_me.append(line.rstrip())

    print "Database: {0}".format(db)
    print "Max Records: {0}".format(retmax)
    print "Output Directory: {0}".format(outdir)
    print "Query: {0}".format(term)
    print "----------"

    print "Searching for records..."
    handle = Entrez.esearch(db=db, retmax=retmax, term=term)
    esearch = Entrez.read(handle)
    print "\tFound {0} records.\n".format(esearch["Count"])

    print "Downloading records..."
    accessions = []
    for uuid in esearch['IdList']:
        handle = Entrez.esummary(db=db, id=uuid)
        esummary = Entrez.read(handle)

        if esummary[0]["Caption"] not in skip_me:
            accessions.append(esummary[0]["Caption"])

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

    print "Outputting list of completed genomes..."
    output = '{0}/completed-genomes.txt'.format(outdir)
    with open(output, 'w') as fh:
        fh.write('\n'.join(accessions))
