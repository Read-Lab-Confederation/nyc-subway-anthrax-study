#! /usr/bin/env python

"""
Author: Robert A Petit III.

Date: 04/04/2016

Read Jellyfish counts of BCG and B. anthracis specific kmersand determine the
BCG kmer coverage and the false positives for B. anthracis kmers.
"""

import argparse as ap
import os
import glob
import subprocess
import numpy as np


def get_total_hits(counts):
    """Return the number of kmers with a hit (not zero)."""
    hits = 0
    for c in counts:
        if int(c):
            hits += 1
    return hits


def trim_zeros(a):
    """Remove zeros from NumPy array, if all zeros set it to 0."""
    return a[a != 0]


def get_stats(counts):
    """Return mean of list of lists."""
    stats = {
        'means': [],
        'medians': [],
        'zero_means': [],
        'zero_medians': []
    }
    for c in counts:
        a = np.array(c)
        no_zero = trim_zeros(a)
        stats['means'].append(np.mean(a))
        stats['medians'].append(np.median(a))
        if len(no_zero):
            stats['zero_means'].append(np.mean(no_zero))
            stats['zero_medians'].append(np.median(no_zero))

    if not len(stats['zero_means']):
        stats['zero_means'] = np.array([0])
        stats['zero_medians'] = np.array([0])

    return stats


def jellyfish_query(sample, kmers, jellyfish, ba=False):
    """Run Jellyfish query against a given sample."""
    cmd = [jellyfish, 'query', '-s', kmers, sample]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    counts = []
    for line in stdout.split('\n'):
        if not line:
            continue
        kmer, count = line.rstrip().split(' ')
        counts.append(int(count))
        if ba and int(count):
            print(sample, kmer)

    return counts


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='get-kmer-coverage.py', conflict_handler='resolve',
        description=("Read Jellyfish counts of non-B. anthracis BCG members "
                     "and determine the BCG kmer coverage and the False "
                     "positive rate for B. anthracis kmers.")
    )

    parser.add_argument('directory', type=str, metavar="JELLYFISH_DIR",
                        help='Directory containing Jellyfish counts.')
    parser.add_argument('output', type=str, metavar="OUTPUT",
                        help='File to write output to.')
    parser.add_argument('bcg', type=str, metavar="BCG_KMERS",
                        help='BCG specific k-mers (FASTA format).')
    parser.add_argument('ba', type=str, metavar="BA_KMERS",
                        help='BA specific k-mers (FASTA format).')
    parser.add_argument('jellyfish', type=str, metavar="JELLYFISH_PATH",
                        help=('Path to Jellyfish executable.'))
    parser.add_argument('--num_cpu', default=1, type=int, metavar="INT",
                        help='Number of processors to use. (Default: 1)')

    args = parser.parse_args()

    results = {}
    # List *.jf file in the given directory
    directory = '{0}/*.jf'.format(args.directory)
    for jf_file in glob.glob(directory):
        sample = os.path.basename(jf_file).replace('.jf', '')
        results[sample] = {'bcg': [], 'ba': []}
        results[sample]['bcg'] = jellyfish_query(jf_file, args.bcg,
                                                 args.jellyfish)
        results[sample]['ba'] = jellyfish_query(jf_file, args.ba,
                                                args.jellyfish, ba=True)

    with open(args.output, 'w') as fh:
        for sample, counts in sorted(results.items()):
            bcg_stats = get_stats(counts['bcg'])
            ba_stats = get_stats(counts['ba'])
            # Sample, BCG kmer hits, BA kmer hits, BCG Count Mean,
            # BA Count Mean, BCG Count Medians, BA Count Medians,
            # BCG Count Mean (no zeros), BA Count Mean (no zeros),
            # BCG Count Median (no zeros), BA Count Median (no zeros)
            fh.write(("{0}\t{1}/{2}\t{3}/{4}\t{5:.6f}\t{6:.6f}\t{7:.6f}\t"
                      "{8:.6f}\t{9:.6f}\t{10:.6f}\t{11:.6f}\t{12:.6f}"
                      "\n").format(
                sample,
                get_total_hits(counts['bcg']), len(counts['bcg']),
                get_total_hits(counts['ba']), len(counts['ba']),
                np.mean(np.array(bcg_stats['means'])),
                np.mean(np.array(ba_stats['means'])),
                np.median(np.array(bcg_stats['medians'])),
                np.median(np.array(ba_stats['medians'])),
                np.mean(np.array(bcg_stats['zero_means'])),
                np.mean(np.array(ba_stats['zero_means'])),
                np.median(np.array(bcg_stats['zero_medians'])),
                np.median(np.array(ba_stats['zero_medians'])),
            ))
