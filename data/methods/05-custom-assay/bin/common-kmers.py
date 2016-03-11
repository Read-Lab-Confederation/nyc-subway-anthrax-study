#! /usr/bin/env python

"""
Author: Robert A Petit III.

Date: 01/06/2016

Read a set of kmer counts created via Jellyfish and output on those kmers
which are found in each member.
"""

import argparse as ap

if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='common-kmers.py', conflict_handler='resolve',
        description="Read Jellyfish kmer counts and output only kmers in all."
    )

    parser.add_argument('input', type=str, metavar="JELLYFISH_DUMPS",
                        help=('List of Jellyfish dumps.'))
    parser.add_argument('rrna', type=str, metavar="RRNA_DUMP",
                        help=('A dump of rRNA kmers.'))
    parser.add_argument('output', type=str, metavar="OUTPUT",
                        help=('File to write common kmers to.'))
    parser.add_argument('--limit', default=0, type=int, metavar="INT",
                        help='Maximum number of files to read (Default: all)')
    parser.add_argument('--single_input', action="store_true",
                        help='The input is not a list on jellyfish counts.')

    args = parser.parse_args()

    common_kmers = None
    first_file = True
    total = 0

    if not args.single_input:
        with open(args.input, 'r') as fh:
            # Read file
            for jf_dump in fh:
                kmers = []
                with open(jf_dump.rstrip(), 'r') as f:
                    for line in f:
                        kmer, count = line.split(' ')
                        kmers.append(kmer)

                total += 1
                if first_file:
                    first_file = False
                    common_kmers = set(kmers)
                else:
                    common_kmers = common_kmers.intersection(set(kmers))

                print("{0}\t{1}".format(total, len(common_kmers)))

                if args.limit > 0 and args.limit == total:
                    break
    else:
        common_kmers = []
        with open(args.input, 'r') as fh:
            for line in fh:
                kmer, count = line.split(' ')
                common_kmers.append(kmer)

    print('Filter rRNA kmers...')
    rrna_kmers = {}
    with open(args.rrna, 'r') as fh:
        for line in fh:
            kmer, count = line.split(' ')
            rrna_kmers[kmer] = True

    total = 0
    with open(args.output, 'w') as fh_out:
        for kmer in sorted(list(common_kmers)):
            if kmer not in rrna_kmers:
                total += 1
                # Write as FASTA, it plays well with jellyfish query
                fh_out.write('>{0}\n{0}\n'.format(kmer))

    print('{0} kmers written'.format(total))
