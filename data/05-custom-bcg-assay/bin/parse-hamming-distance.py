#! /usr/bin/env python

"""
Author: Robert A Petit III.

Date: 03/28/2016

Read Jellyfish counts of non-B. anthracis BCG members and determine the BCG
kmer coverage and the false positives for B. anthracis kmers.
"""

import argparse as ap
import json

if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='parse-hamming-distance.py', conflict_handler='resolve',
        description=("Read BLAST JSON output and determine the hamming "
                     "distance between a k-mer and its best hit.")
    )

    parser.add_argument('json', type=str, metavar="BLAST_JSON",
                        help=('BLAST output in JSON format.'))
    args = parser.parse_args()

    json_data = []
    with open(args.json) as fh:
        record = []
        first_line = True
        for line in fh:
            if line.startswith('{') and not first_line:
                json_data.append(json.loads(''.join(record)))
                record = []
            if first_line:
                first_line = False
            record.append(line)
        json_data.append(json.loads(''.join(record)))

    for entry in json_data:
        hit = entry['BlastOutput2']['report']['results']['search']
        hsp = hit['hits'][0]['hsps'][0]

        # Includes mismatches and gaps
        mismatch = hsp['align_len'] - hsp['identity']

        # Hamming distance
        hd = mismatch
        if hit['query_len'] > hsp['align_len']:
            # Include those bases that weren't aligned
            hd = hit['query_len'] - hsp['align_len'] + mismatch

        print('{0}\t{1}'.format(
            hit['query_title'],
            hd
        ))
