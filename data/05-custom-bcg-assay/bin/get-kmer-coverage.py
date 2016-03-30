#! /usr/bin/env python

"""
Author: Robert A Petit III.

Date: 03/28/2016

Read a kmer-by-sample files and determine the per sample and per kmer
coverages.
"""

import argparse as ap
import os
import glob
import numpy as np


def gziplines(fname):
    """Use zcat to parse gzip file."""
    from subprocess import Popen, PIPE
    f = Popen(['zcat', fname], stdout=PIPE)
    for line in f.stdout:
        yield line

if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='get-kmer-coverage.py', conflict_handler='resolve',
        description=("Read a kmer-by-sample files and determine the per "
                     "sample and per kmer coverages.")
    )

    parser.add_argument('directory', type=str, metavar="SIMUALTION_DIR",
                        help=('Directory containing coverage simulations with '
                              '*-kmer-by-sample.txt.gz files to processed.'))
    parser.add_argument('samples', type=str, metavar="SAMPLES",
                        help='Samples with organism name')
    parser.add_argument('filter', type=str, metavar="FILTER",
                        help='Organism name to filter results by.')
    parser.add_argument('per_sample', type=str, metavar="PER_SAMPLE_STATS",
                        help='File to write per sample stats to.')
    parser.add_argument('per_kmer', type=str, metavar="PER_KMER_STATS",
                        help='File to write per kmer stats to.')

    args = parser.parse_args()

    organism_name = {}
    with open(args.samples, 'r') as fh:
        for line in fh:
            sample, organism = line.rstrip().split('\t')
            organism_name[sample] = organism.lower()

    per_sample = {}
    per_kmer = {}
    for coverage in os.walk(args.directory).next()[1]:
        directory = '{0}/{1}/*-kmer-by-sample.txt.gz'.format(
            args.directory, coverage
        )
        # Expected coverage
        cov = coverage.replace('x', '')
        per_sample[cov] = {}
        per_kmer[cov] = {}
        for input_file in glob.glob(directory):
            for line in gziplines(input_file):
                sample, kmer, count = line.rstrip().split('\t')
                if sample not in per_sample[cov]:
                    per_sample[cov][sample] = [int(count)]
                else:
                    per_sample[cov][sample].append(int(count))

                if kmer not in per_kmer[cov]:
                    per_kmer[cov][kmer] = {
                        'filter': [],
                        'non-filter': []
                    }

                if args.filter.lower() in organism_name[sample]:
                    per_kmer[cov][kmer]['filter'].append(int(count))
                else:
                    per_kmer[cov][kmer]['non-filter'].append(int(count))

    with open(args.per_sample, 'w') as fh:
        fh.write('\t'.join([
            'sample', 'filter', 'expected_coverage', 'mean', 'median',
            'mean_nonzero', 'median_nonzero'
        ]))
        fh.write('\n')
        for coverage, samples in sorted(per_sample.items()):
            for sample, counts in samples.items():
                a = np.array(counts)
                a_nonzero = a[a != 0]
                if not len(a_nonzero):
                    a_nonzero = np.array([0])
                fh.write('{0}\t{1}\t{2}\t{3:.4f}\t{4}\t{5:.4f}\t{6}\n'.format(
                    sample,
                    args.filter.lower() in organism_name[sample],
                    coverage,
                    np.mean(a),
                    np.median(a),
                    np.mean(a_nonzero),
                    np.median(a_nonzero)
                ))

    with open(args.per_kmer, 'w') as fh:
        fh.write('\t'.join([
            'kmer', 'expected_coverage', 'filtered_mean', 'non_filtered_mean',
            'filtered_median', 'non_filtered_median', 'filtered_mean_nonzero',
            'non_filtered_mean_nonzero', 'filtered_median_nonzero',
            'non_filtered_median_nonzero'
        ]))
        fh.write('\n')
        for coverage, kmers in sorted(per_kmer.items()):
            for kmer, counts in sorted(kmers.items()):
                if not len(counts['filter']):
                    counts['filter'].append(0)
                filtered = np.array(counts['filter'])
                filtered_nonzero = filtered[filtered != 0]
                if not len(filtered_nonzero):
                    filtered_nonzero = np.array([0])

                if not len(counts['non-filter']):
                    counts['non-filter'].append(0)
                non_filtered = np.array(counts['non-filter'])

                non_filtered_nonzero = non_filtered[non_filtered != 0]
                if not len(non_filtered_nonzero):
                    non_filtered_nonzero = np.array([0])

                fh.write(('{0}\t{1}\t{2:.4f}\t{3:.4f}\t{4}\t{5}\t{6:.4f}\t'
                          '{7:.4f}\t{8}\t{9}\n').format(
                    kmer,
                    coverage,
                    np.mean(filtered),
                    np.mean(non_filtered),
                    np.median(filtered),
                    np.median(non_filtered),
                    np.mean(filtered_nonzero),
                    np.mean(non_filtered_nonzero),
                    np.median(filtered_nonzero),
                    np.median(non_filtered_nonzero)
                ))
