#! /usr/bin/env python
""" Organism counts of top blast hits from reads aligned to toxin genes. """
import fnmatch
import operator
import os
import re

hits = []
for root, dirnames, filenames in os.walk('data/results'):
    for filename in fnmatch.filter(filenames, 'all.summary'):
        hits.append(os.path.join(root, filename))

pattern = re.compile(r"data/results/02-mapping/(?P<study>.*)/(?P<sample>.*)/"
                     "(?P<plasmid>.*)/aligned-genes/all.summary")
results = {}
for hit in hits:
    m = pattern.search(hit)
    study = m.group('study')
    plasmid = m.group('plasmid')
    sample = m.group('sample')

    if study not in results:
        results[study] = {}

    for line in open(hit, "r"):
        line = line.strip()
        count, organism = line.split(" ", 1)

        if organism not in results[study]:
            results[study][organism] = int(count)
        else:
            results[study][organism] += int(count)


for study in results:
    output = "data/results/02-mapping/{0}/top_blast_hit_counts.txt".format(
        study
    )
    out = open(output, "w")

    sorted_vals = sorted(results[study].items(), key=operator.itemgetter(1),
                         reverse=True)
    for (organism, count) in sorted_vals:
        out.write("{0}\t{1}\n".format(organism, count))

    out.write("\n")
    out.close()
