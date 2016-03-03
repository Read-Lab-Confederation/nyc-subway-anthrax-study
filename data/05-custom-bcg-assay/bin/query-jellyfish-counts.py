#! /usr/bin/env python
"""Count kmers."""
from __future__ import print_function
from os.path import basename, splitext
import argparse as ap
import glob
import sys
import subprocess
import numpy as np


def jellyfish_query(sample, kmers, jellyfish):
    """Run Jellyfish query against a given sample."""
    cmd = [jellyfish, 'query', '-s', kmers, sample]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return stdout.split('\n')


def print_row(string, out=sys.stdout):
    """Print the row values to stdout, or stderr if specified."""
    print(string, file=out)

if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='download-genomes.py', conflict_handler='resolve',
        description="Download sequences individually."
    )

    parser.add_argument('kmers', type=str, metavar="KMERS",
                        help=('Kmers (FASTA format) to search for.'))
    parser.add_argument('jellyfish_dir', type=str, metavar="JELLYFISH_DIR",
                        help=('Directory of Jellyfish counts.'))
    parser.add_argument('jellyfish', type=str, metavar="JELLYFISH_PATH",
                        help=('Path to Jellyfish executable.'))
    parser.add_argument('output_stats', type=str, metavar="OUTPUT_STATS",
                        help=('File to output per sample stats.'))
    parser.add_argument('output_kmer', type=str, metavar="OUTPUT_KMER",
                        help=('File to output kmers per sample.'))
    parser.add_argument('output_kmers', type=str, metavar="OUTPUT_KMERS",
                        help=('File to output kmer counts.'))

    args = parser.parse_args()
    kmers = []
    stats = []
    kmer_by_sample = {}
    kmer_counts = {}

    cols = ('Sample\tTotal Queried Kmers\tHits\tPercent\tSingleton Hits\t'
            'Singleton Percent\tHits (excluding singletons)\t'
            'Percent(excluding singletons)\tMin Count\tMean Count\t'
            'Median Count\tMax Count')
    stats.append(cols)
    print_row(cols, sys.stderr)

    directory = '{0}/*.jf'.format(args.jellyfish_dir)
    for input_file in glob.glob(directory):

        sample = splitext(basename(input_file))[0]
        kmer_by_sample[sample] = {}
        counts = []
        hit = 0
        singletons = 0

        for line in jellyfish_query(input_file, args.kmers, args.jellyfish):
            if not line:
                continue

            try:
                kmer, count = line.rstrip().split(' ')
                count = int(count)
                counts.append(count)
            except ValueError:
                print('LINE: {0}'.format(line))
                print('COUNTS: {0}'.format(len(counts)))
                sys.exit()

            if count:
                kmer_by_sample[sample][kmer] = count
                if kmer not in kmer_counts:
                    kmer_counts[kmer] = 1
                else:
                    kmer_counts[kmer] += 1

                if count == 1:
                    singletons += 1
                hit += 1

        a = np.array(counts)

        row = ('{0}\t{1}\t{2}\t{3:.2f}\t{4}\t{5:.2f}\t{6}\t{7:.2f}\t{8}\t'
               '{9:.2f}\t{10:.2f}\t{11}').format(
            sample,
            len(counts),
            hit,
            float(hit) / len(counts),
            singletons,
            float(singletons) / hit if (hit) else 0.00,
            hit - singletons,
            float(hit - singletons) / len(counts),
            min(counts),
            np.mean(a),
            np.median(a),
            max(counts)
        )
        stats.append(row)

        print_row(row, sys.stderr)

    with open(args.output_stats, 'w') as fh:
        for row in stats:
            fh.write("{0}\n".format(row))

    with open(args.output_kmer, 'w') as fh:
        for sample in kmer_by_sample:
            for kmer in kmer_by_sample[sample]:
                fh.write('{0}\t{1}\t{2}\n'.format(
                    sample,
                    kmer,
                    kmer_by_sample[sample][kmer]
                ))

    with open(args.output_kmers, 'w') as fh:
        for kmer, count in kmer_counts.items():
            fh.write('{0}\t{1}\n'.format(
                kmer, count
            ))
