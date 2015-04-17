#! /usr/bin/env python
""" Summary statistics of plasmid coverage. """
import fnmatch
import gzip
import os
import re

import numpy as np

coverages = []
for root, dirnames, filenames in os.walk('data/results'):
    for filename in fnmatch.filter(filenames, '*.coverage.gz'):
        coverages.append(os.path.join(root, filename))

pattern = re.compile(r"data/results/02-mapping/(?P<study>.*)/(?P<sample>.*)/"
                     "(?P<plasmid>.*)/coverage/(?P<coverage>.*).coverage.gz")

results = {}
for coverage in coverages:
    m = pattern.search(coverage)
    study = m.group('study')
    plasmid = m.group('plasmid')

    if study not in results:
        results[study] = {}

    if plasmid not in results[study]:
        results[study][plasmid] = []

    fh = gzip.open(coverage, 'rb')
    per_base_coverage = np.array([int(l.split('\t')[2].rstrip()) for l in fh])
    if len(per_base_coverage) > 0:
        results[study][plasmid].append([
            m.group('sample'),
            np.min(per_base_coverage),
            np.percentile(per_base_coverage, 25),
            np.median(per_base_coverage),
            np.mean(per_base_coverage),
            np.percentile(per_base_coverage, 75),
            np.max(per_base_coverage)
        ])
    else:
        results[study][plasmid].append([m.group('sample'), 0, 0, 0, 0, 0, 0])

for study in results:
    output = "data/results/02-mapping/{0}/coverage_summary.txt".format(study)
    out = open(output, "w")
    for plasmid in results[study]:
        out.write("{0}\n".format(plasmid))
        out.write("Study\tMin.\t1st Qu.\tMedian\tMean\t3rd Qu.\tMax\n")
        for i in results[study][plasmid]:
            row = '{0}\t{1}\t{2:.4f}\t{3:.4f}\t{4:.4f}\t{5:.4f}\t{6}'.format(
                i[0], i[1], i[2], i[3], i[4], i[5], i[6]
            )
            out.write("{0}\n".format(row))
        out.write("\n\n")
    out.close()
