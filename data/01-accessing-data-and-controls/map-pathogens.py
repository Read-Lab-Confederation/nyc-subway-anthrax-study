#! /usr/bin/env python
"""
Create symbolic links of yersinia and anthracis in a new folder.

sys.argv[1]: data/pathogens-to-samples.txt
sys.argv[2]: data/runs-to-samples.txt

mkdir sra-pathogens
mkdir sra-pathogens/anthracis
mkdir sra-pathogens/yersinia
"""
import os
import errno
import glob
import sys


def mkdir(directory):
    """ Make a directory using sample id and pthogen. """
    if not os.path.exists(directory):
        try:
            print 'Making {0}'.format(directory)
            os.makedirs(directory)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def symbolic_link(indir, outdir, run_acc, sample_id):
    """ Create symlink of the SRA fastq, and rename it with sample ID. """
    for run_file in glob.glob(indir):
        run_file = os.path.abspath(run_file)
        sample_file = outdir + os.path.basename(run_file)
        print 'Creating symlink: {0} -> {1}'.format(run_file, sample_file)
        os.symlink(run_file, sample_file)


# Make some directories
mkdir('sra-pathogens')
mkdir('sra-pathogens/anthracis')
mkdir('sra-pathogens/yersinia')

# Read sample ids and the pathogen
sample_ids = {}
fh = open(sys.argv[1], 'r')
for line in fh:
    line = line.rstrip()
    cols = line.split('\t')

    # cols[0]: Sample ID, cols[1]: Pathogen (anthracis|yersinia)
    sample_ids[cols[0]] = cols[1]
fh.close()

# Read sample ids and their run accessions
runs_to_samples = {}
fh = open(sys.argv[2], 'r')
for line in fh:
    line = line.rstrip()
    cols = line.split('\t')

    # cols[0]: Run accession, cols[1]: Sample ID
    if cols[1] in sample_ids:
        indir = 'sra-fastq/{0}/*'.format(cols[0])
        outdir = 'sra-pathogens/{0}/{1}/'.format(sample_ids[cols[1]], cols[1])
        mkdir(outdir)
        symbolic_link(indir, outdir, cols[0], cols[1])
fh.close()
