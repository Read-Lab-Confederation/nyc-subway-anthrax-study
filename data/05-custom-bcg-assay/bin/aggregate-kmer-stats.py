#! /usr/bin/env python
"""Get column means grouped by kmer."""
import argparse as ap
import pandas as pd
import numpy as np


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='aggregate-kmer-stats.py', conflict_handler='resolve',
        description="Assess the results of all the simulations."
    )

    parser.add_argument('results', type=str, metavar="VALIDATION_OUTPUT",
                        help='Output the mean kmer stats of all simulations.')
    args = parser.parse_args()

    df = pd.read_csv(args.results, sep='\t')
    df.columns = ['Kmer', 'Total Samples', 'Expected TP', 'Expected TN',
                  'TP', 'FN', 'FP', 'TN', 'Accuracy', 'Precision',
                  'Specificity', 'Sensitivity', 'FPR']

    grouped = df.groupby(['Kmer'], as_index=False)
    means = grouped.aggregate(np.mean)

    output = args.results.replace('.txt', '-aggregated.txt')
    means.to_csv(output, sep='\t')
