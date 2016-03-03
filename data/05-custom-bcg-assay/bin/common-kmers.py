#! /usr/bin/env python

"""
Author: Robert A Petit III.

Date: 01/06/2016

Read a set of kmer counts created via Jellyfish and output on those kmers
which are found in each member.
"""
import argparse as ap
import glob
import subprocess
import tempfile


def jellyfish_dump(sample, jellyfish):
    """Run Jellyfish query against a given sample."""
    cmd = [jellyfish, 'dump', '-c', sample]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return stdout.split('\n')


def jellyfish_query(sample, kmers, jellyfish):
    """Run Jellyfish query against a given sample."""
    stdout = None
    try:
        tmp = tempfile.NamedTemporaryFile()
        for kmer in kmers:
            tmp.write('>{0}\n{0}\n'.format(kmer))

        p = subprocess.Popen(
            [jellyfish, 'query', '-s', tmp.name, sample],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = p.communicate()

    finally:
        tmp.close()

    return stdout.split('\n')


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='common-kmers.py', conflict_handler='resolve',
        description="Read Jellyfish kmer counts and output only kmers in all."
    )

    parser.add_argument('jellyfish_dir', type=str, metavar="JELLYFISH_DIR",
                        help=('Directory of Jellyfish counts.'))
    parser.add_argument('jellyfish', type=str, metavar="JELLYFISH_PATH",
                        help=('Path to Jellyfish executable.'))
    parser.add_argument('rrna', type=str, metavar="RRNA_JF",
                        help=('A dump of rRNA kmers.'))
    parser.add_argument('output', type=str, metavar="OUTPUT",
                        help=('File to write common kmers to.'))
    parser.add_argument('--limit', default=0, type=int, metavar="INT",
                        help='Maximum number of files to read (Default: all)')
    parser.add_argument('--single_input', action="store_true",
                        help='The input is not a list on jellyfish counts.')

    args = parser.parse_args()

    common_kmers = None
    first_file = True
    total = 0

    if not args.single_input:
        directory = '{0}/*.jf'.format(args.jellyfish_dir)
        for input_file in glob.glob(directory):
            total += 1
            kmers = []
            if first_file:
                first_file = False
                for line in jellyfish_dump(input_file, args.jellyfish):
                    if not line:
                        continue

                    kmer, count = line.split(' ')
                    kmers.append(kmer)
                common_kmers = set(kmers)
            else:
                for line in jellyfish_query(input_file, list(common_kmers),
                                            args.jellyfish):
                    if not line:
                        continue

                    kmer, count = line.rstrip().split(' ')

                    if int(count):
                        kmers.append(kmer)

                common_kmers = common_kmers.intersection(set(kmers))

            print("{0}\t{1}".format(total, len(common_kmers)))

            if args.limit > 0 and args.limit == total:
                break
    else:
        common_kmers = []
        with open(args.jellyfish_dir, 'r') as fh:
            for line in fh:
                kmer, count = line.split(' ')
                common_kmers.append(kmer)

    print('Filter rRNA kmers...')
    rrna_kmers = {}
    for line in jellyfish_dump(args.rrna, args.jellyfish):
        if not line:
            continue
        kmer, count = line.split(' ')
        rrna_kmers[kmer] = True

    total = 0
    with open(args.output, 'w') as fh_out:
        for kmer in sorted(list(common_kmers)):
            if kmer not in rrna_kmers:
                total += 1
                # Write as FASTA, it plays well with jellyfish query
                fh_out.write('>{0}\n{0}\n'.format(kmer))

    print('{0} kmers written'.format(total))
