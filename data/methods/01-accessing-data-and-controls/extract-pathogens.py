#! /usr/bin/env python
"""
Read supplementary table and extract pathogens.

sys.argv[1]: data/DataTable5-metaphlan-metadata_v19.txt
Extract the sample id and the columns which pertain to yersinia and
anthracis.
"""


def split_header(header):
    """
    Some headers are really long, return only the last portion.

    Example Headers:
    k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillaceae|g__Bacillus|s__Bacillus_anthracis
    k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillaceae|g__Bacillus|s__Bacillus_anthracis|t__Bacillus_anthracis_unclassified
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Yersinia
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Yersinia|s__Yersinia_unclassified

    Return Headers:
    s__Bacillus_anthracis
    t__Bacillus_anthracis_unclassified
    g__Yersinia
    s__Yersinia_unclassified
    """
    return header.split('|')[-1]


import sys
fh = open(sys.argv[1], 'rU')
pathogens = []
pathogen_name = {}
total = 0
for line in fh:
    line = line.rstrip()
    cols = line.split('\t')
    if len(pathogens) == 0:
        for i in xrange(len(cols)):
            if "anthracis" in cols[i] or "Yersinia" in cols[i]:
                pathogens.append(i)
                pathogen_name[i] = cols[i]

        print '\t'.join([
            cols[0],
            'pathogen',
            split_header(cols[pathogens[0]]),
            split_header(cols[pathogens[1]]),
            split_header(cols[pathogens[2]]),
            split_header(cols[pathogens[3]])
        ])
    else:
        sample_contains_pathogen = False
        is_anthracis = False
        is_yersinia = False
        for i in pathogens:
            if float(cols[i]) > 0:
                if "anthracis" in pathogen_name[i]:
                    is_anthracis = True
                elif "Yersinia" in pathogen_name[i]:
                    is_yersinia = True
                sample_contains_pathogen = True

        if sample_contains_pathogen:
            pathogen = None
            if is_anthracis and is_yersinia:
                pathogen = 'anthracis/yersinia'
            elif is_anthracis:
                pathogen = 'anthracis'
            elif is_yersinia:
                pathogen = 'yersinia'

            print '\t'.join([cols[0], pathogen, cols[pathogens[0]],
                            cols[pathogens[1]], cols[pathogens[2]],
                            cols[pathogens[3]]])
fh.close()
