#! /usr/bin/env python

"""
Author: Robert A Petit III.

Date: 03/28/2016

Read Jellyfish counts of non-B. anthracis BCG members and determine the BCG
kmer coverage and the false positives for B. anthracis kmers.
"""

import argparse as ap
import os
import subprocess
import multiprocessing as mp
import numpy as np


def trim_zeros(a):
    """Remove zeros from NumPy array, if all zeros set it to 0."""
    a_nonzero = a[a != 0]
    if not len(a_nonzero):
        a_nonzero = np.array([0])

    return a_nonzero


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
        stats['zero_means'].append(np.mean(no_zero))
        stats['zero_medians'].append(np.median(no_zero))

    return stats


def query(files, kmers, jf, num_cpu):
    """Multiprocess the Jellyfish queries."""
    pool = mp.Pool(processes=num_cpu)
    results = [pool.apply_async(jellyfish_query, args=(f, kmers, jf))
               for f in files]
    pool.close()
    pool.join()

    return [p.get() for p in results]


def jellyfish_query(sample, kmers, jellyfish):
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

    return counts


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='get-kmer-coverage.py', conflict_handler='resolve',
        description=("Read Jellyfish counts of non-B. anthracis BCG members "
                     "and determine the BCG kmer coverage and the False "
                     "positive rate for B. anthracis kmers.")
    )

    parser.add_argument('directory', type=str, metavar="SIMUALTION_DIR",
                        help=('Directory containing coverage simulations with '
                              'Jellyfish counts.'))
    parser.add_argument('samples', type=str, metavar="SAMPLES",
                        help='Samples with organism name')
    parser.add_argument('bcg', type=str, metavar="BCG_KMERS",
                        help='BCG specific k-mers (FASTA format).')
    parser.add_argument('ba', type=str, metavar="BA_KMERS",
                        help='BA specific k-mers (FASTA format).')
    parser.add_argument('jellyfish', type=str, metavar="JELLYFISH_PATH",
                        help=('Path to Jellyfish executable.'))
    parser.add_argument('--num_cpu', default=1, type=int, metavar="INT",
                        help='Number of processors to use. (Default: 1)')

    args = parser.parse_args()

    # Read sample file and filter out B. anthracis samples
    accessions = {}
    with open(args.samples, 'r') as fh:
        for line in fh:
            sample, organism = line.rstrip().split('\t')
            if 'anthracis' not in organism.lower():
                accessions[sample] = organism.lower()

    results = {}
    # Traverse the simulation directory
    for coverage in os.walk(args.directory).next()[1]:
        # Expected coverage
        cov = coverage.replace('x', '')
        results[cov] = {}
        for sample, organism in sorted(accessions.items()):
            if sample not in results[cov]:
                results[cov][sample] = {'bcg': [], 'ba': []}
            # Traverse simulation runs
            files = []
            sim_dir = '{0}/{1}/'.format(args.directory, coverage)
            for sim in os.walk(sim_dir).next()[1]:
                files.append('{0}/{1}/{2}/{3}.jf'.format(
                    args.directory, coverage, sim, sample
                ))
            results[cov][sample]['bcg'] = query(files, args.bcg,
                                                args.jellyfish, args.num_cpu)
            results[cov][sample]['ba'] = query(files, args.ba,
                                               args.jellyfish, args.num_cpu)

    for cov, samples in sorted(results.items()):
        for sample, counts in sorted(samples.items()):
            bcg_stats = get_stats(counts['bcg'])
            ba_stats = get_stats(counts['ba'])

            # Coverage, Sample, BCG Count Median, BA Count Median,
            # BCG Count Median (no zeros), BA Count Median (no zeros)
            print(("{0}\t{1}\t{2:.6f}\t{3:.6f}\t{4:.6f}\t{5:.6f}\t{6:.6f}\t"
                   "{7:.6f}\t{8:.6f}\t{9:.6f}").format(
                cov,
                sample,
                np.mean(np.array(bcg_stats['means'])),
                np.mean(np.array(ba_stats['means'])),
                np.median(np.array(bcg_stats['medians'])),
                np.median(np.array(ba_stats['medians'])),
                np.mean(np.array(bcg_stats['zero_means'])),
                np.mean(np.array(ba_stats['zero_means'])),
                np.median(np.array(bcg_stats['zero_medians'])),
                np.median(np.array(ba_stats['zero_medians'])),
            ))
