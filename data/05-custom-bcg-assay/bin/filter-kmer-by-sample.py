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
        prog='download-genomes.py', conflict_handler='resolve',
        description="Download sequences individually."
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
        sample, kmer = line.split('\t')[0:2]
        sample = sample.split('-')[0]
        if kmer not in kmers:
            kmers[kmer] = {
                'TP': [],
                'FP': []
            }

        if args.filter.lower() in organism_name[sample]:
            kmers[kmer]['TP'].append(sample)
        else:
            kmers[kmer]['FP'].append(sample)

    for kmer in sorted(kmers):
        tp = len(set(kmers[kmer]['TP']))
        fn = total_tp - tp
        fp = len(set(kmers[kmer]['FP']))
        tn = total_tn - fp

        # See Olson et al. 2015 frontiers in Genetics for calculations
        accuracy = float(tp + tn) / (tp + fp + fn + tn)
        specificity = float(tn) / (fp + tn)
        sensitivity = float(tp) / (tp + fn)
        precision = float(tp) / (tp + fp)
        fpr = float(fp) / (tn + fp)

        print(('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8:.2f}\t{9:.2f}\t'
               '{10:.2f}\t{11:.2f}\t{12:.2f}').format(
            kmer,
            total_tp + total_tn,
            total_tp,
            total_tn,
            tp,
            fn,
            fp,
            tn,
            accuracy,
            specificity,
            sensitivity,
            precision,
            fpr
        ))
