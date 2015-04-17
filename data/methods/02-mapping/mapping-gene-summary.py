#! /usr/bin/env python
""" Total number of reads aligned to toxin related genes. """
import fnmatch
import os
import re

genes = []
for root, dirnames, filenames in os.walk('data/results'):
    for filename in fnmatch.filter(filenames, '*.fasta'):
        genes.append(os.path.join(root, filename))

pattern = re.compile(r"data/results/02-mapping/(?P<study>.*)/(?P<sample>.*)/"
                     "(?P<plasmid>.*)/aligned-genes/(?P<gene>.*).fasta")

results = {}
for gene in genes:
    print gene
    m = pattern.search(gene)
    study = m.group('study')
    plasmid = m.group('plasmid')
    sample = m.group('sample')

    if study not in results:
        results[study] = {}

    if plasmid not in results[study]:
        results[study][plasmid] = {
            'genes': []
        }

    if sample not in results[study][plasmid]:
        results[study][plasmid][sample] = {}

    read_count = 0
    for line in open(gene, "r"):
        if line.startswith(">"):
            read_count += 1

    if m.group('gene') not in results[study][plasmid]['genes']:
        results[study][plasmid]['genes'].append(m.group('gene'))
    results[study][plasmid][sample][m.group('gene')] = str(read_count)

for study in results:
    output = "data/results/02-mapping/{0}/toxin_read_counts.txt".format(study)
    out = open(output, "w")
    for plasmid in results[study]:
        out.write("{0}\n".format(plasmid))
        out.write("Sample\t{0}\n".format(
            '\t'.join(sorted(results[study][plasmid]['genes']))
        ))
        for sample in results[study][plasmid]:
            if sample != "genes":
                read_counts = [sample]
                for gene in sorted(results[study][plasmid]['genes']):
                    read_counts.append(results[study][plasmid][sample][gene])

                out.write("{0}\n".format('\t'.join(read_counts)))
        out.write("\n")
    out.close()
