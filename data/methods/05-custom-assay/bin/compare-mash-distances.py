#! /usr/bin/env python

"""
Author: Robert A Petit III.

Date: 02/24/2016

Determine if a genome should be include in a given group. This essentially uses
mash to determine if a nother Bacillus genome should be included within the
Bacillus cereus Group. If a genome's maximum distance against all members of
the group less than the maximum self hit distances it should be included in
the group.
"""

import argparse as ap
from os.path import basename, splitext

if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='compare-mash-distances.py', conflict_handler='resolve',
        description="Use Mash distances to determine group inclusion."
    )

    parser.add_argument('group', type=str, metavar="MASH_GROUP",
                        help='Mash distances of the group compared to itself.')
    parser.add_argument('mash_list', type=str, metavar="MASH_LIST",
                        help=('List of mash distances compared to the group '
                              'of interest.'))
    parser.add_argument('output', type=str, metavar="OUTPUT",
                        help=('File to write new members to.'))

    args = parser.parse_args()

    group_distances = []
    with open(args.group, 'r') as fh:
        for line in fh:
            group_distances.append(float(line.split('\t')[2]))

    with open(args.output, 'w') as fh_out, open(args.mash_list, 'r') as fh:
        # Read file
        for mash_file in fh:
            distances = []
            with open(mash_file.rstrip(), 'r') as f:
                for line in f:
                    distances.append(float(line.split('\t')[2]))

            if max(distances) <= max(group_distances):
                # Genome should be included in the group
                genome = splitext(basename(mash_file))[0]
                print("{0} should be in the group ({1} <= {2})".format(
                    genome, max(distances), max(group_distances)
                ))
                fh_out.write("{0}\n".format(genome))
