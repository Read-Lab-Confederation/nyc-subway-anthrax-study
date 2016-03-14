#! /usr/bin/env python

"""
Author: Robert A Petit III.

Date: 02/26/2016

Read a kmer-by-sample file and filter based on organism.
"""

import argparse as ap


def gziplines(fname):
    """Use zcat to parse gzip file."""
    from subprocess import Popen, PIPE
    f = Popen(['zcat', fname], stdout=PIPE)
    for line in f.stdout:
        yield line

if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='filter-kmer-by-sample.py', conflict_handler='resolve',
        description="Read a kmer-by-sample file and filter based on organism.."
    )

    parser.add_argument('kmers', type=str, metavar="KMER_BY_SAMPLE_GZIP",
                        help=('Gzipped text file of kmer by sample.'))
    parser.add_argument('samples', type=str, metavar="SAMPLES",
                        help=('Samples with organism name'))
    parser.add_argument('filter', type=str, metavar="FILTER",
                        help=('Organism name to filter results by.'))

    args = parser.parse_args()

    organism_name = {}
    total_tp = 0
    total_tn = 0
    with open(args.samples, 'r') as fh:
        for line in fh:
            sample, organism = line.rstrip().split('\t')
            organism_name[sample] = organism.lower()

            if args.filter.lower() in organism_name[sample]:
                total_tp += 1
            else:
                total_tn += 1

    kmers = {}
    for line in gziplines(args.kmers):
        sample, kmer, count = line.split('\t')
        sample = sample.split('-')[0]
        if kmer not in kmers:
            kmers[kmer] = {
                'TP': [],
                'FP': [],
                'FN': [],
                'TN': []
            }

        if args.filter.lower() in organism_name[sample]:
            if int(count):
                kmers[kmer]['TP'].append(sample)
            else:
                kmers[kmer]['FN'].append(sample)
        else:
            if int(count):
                kmers[kmer]['FP'].append(sample)
            else:
                kmers[kmer]['TN'].append(sample)

    for kmer in sorted(kmers):
        tp = len(set(kmers[kmer]['TP']))
        fn = len(set(kmers[kmer]['FN']))
        fp = len(set(kmers[kmer]['FP']))
        tn = len(set(kmers[kmer]['TN']))

        # See Olson et al. 2015 frontiers in Genetics for calculations
        accuracy = float(tp + tn) / (tp + fp + fn + tn)
        specificity = float(tn) / (fp + tn)
        sensitivity = float(tp) / (tp + fn)
        precision = float(tp) / (tp + fp) if tp else 0.00
        fpr = float(fp) / (tn + fp)

        print(('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8:.4f}\t{9:.4f}\t'
               '{10:.4f}\t{11:.4f}\t{12:.4f}').format(
            kmer,
            total_tp + total_tn,
            total_tp,
            total_tn,
            tp,
            fn,
            fp,
            tn,
            accuracy,
            precision,
            specificity,
            sensitivity,
            fpr
        ))
