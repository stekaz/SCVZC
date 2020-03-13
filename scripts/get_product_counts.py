#!/usr/bin/env python3
# coding: utf-8
###############################################################################
#
#    get_product_counts.py
#
#    Counts alignments for products listed in a GenBank file
#
#    Copyright (C) 2020 QIMR Berghofer Medical Research Institute
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

import argparse
import csv
import logging
import os
import sys
import pysam

from contextlib import ExitStack

from Bio import SeqIO

logger = logging.getLogger()


def get_option_parser():

    def formatter_class(prog):
        kwargs = { 'max_help_position': 80, 'width': 100 }
        return argparse.RawTextHelpFormatter(prog, **kwargs)

    parser = argparse.ArgumentParser(formatter_class=formatter_class)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        'gb_file',
        help='Path to the input GenBank file',
        metavar='<GENBANK>')
    group.add_argument(
        'bam_files',
        help='Path to the input BAM file',
        metavar='<BAM>',
        nargs='+')

    group = parser.add_argument_group('output options')
    group.add_argument(
        '-o',
        help='Write output to FILE',
        dest='outfile',
        default='-',
        metavar='FILE')

    group = parser.add_argument_group('filter options')
    group.add_argument(
        '-q',
        help='Only include reads with mapping quality >= INT [0]',
        dest='mapping_quality',
        metavar='INT',
        type=int)
    group.add_argument(
        '-f',
        help='Only include reads with all of the FLAGs in INT present [0]',
        dest='include_flags',
        default=[],
        action='append',
        metavar='INT',
        type=int)
    group.add_argument(
        '-F',
        help='Only include reads with none of the FLAGS in INT present [0]',
        dest='exclude_flags',
        default=[],
        action='append',
        metavar='INT',
        type=int)

    group = parser.add_argument_group('logging options')
    group.add_argument(
        '--debug',
        help='Set logging level to DEBUG',
        action="store_const",
        dest="loglevel",
        const=logging.DEBUG,
        default=logging.WARNING,
    )
    group.add_argument(
        '--verbose',
        help='Set logging level to INFO',
        action="store_const",
        dest="loglevel",
        const=logging.INFO,
    )

    return parser


class MappingQualityFilter:
    def __init__(self, quality):
        self.quality = quality

    def __call__(self, alignment):
        return alignment.mapping_quality >= self.quality


class IncludeFlagFilter:
    def __init__(self, flag):
        self.flag = flag

    def __call__(self, alignment):
        return (alignment.flag & self.flag) == self.flag


class ExcludeFlagFilter:
    def __init__(self, flag):
        self.flag = flag

    def __call__(self, alignment):
        return (alignment.flag & self.flag) == 0


def main():

    parser = get_option_parser()
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s: %(message)s', level=args.loglevel)

    samples = args.bam_files

    # try to remove a common dirname
    if len(set(os.path.split(bam)[0] for bam in args.bam_files)) == 1:
        samples = [os.path.split(bam)[1] for bam in args.bam_files]
        # try to remove a common file extension
        if len(set(sample.split('.', 1)[1] for sample in samples)) == 1:
            samples = [sample.split('.', 1)[0] for sample in samples]

    output_fields = [ "region", "strand", "type", "product" ]

    filters = []

    if args.mapping_quality:
        logging.info("Applying MAPQ filter: {}".format(str(args.mapping_quality)))
        filters.append(MappingQualityFilter(args.mapping_quality))

    for flag in args.include_flags:
        logging.info("Only including reads with all of the FLAGs in: {}".format(str(flag)))
        filters.append(IncludeFlagFilter(flag))

    for flag in args.exclude_flags:
        logging.info("Only including reads with none of the FLAGs in: {}".format(str(flag)))
        filters.append(ExcludeFlagFilter(flag))

    read_callback = lambda aln: all(fltr(aln) for fltr in filters)

    with ExitStack() as stack:

        logging.debug("Opening GenBank file: {}".format(args.gb_file))
        gb_fh = stack.enter_context(open(args.gb_file, 'rU'))

        bams = []

        for bam in args.bam_files:
            logging.debug("Opening BAM file: {}".format(bam))
            bams.append(stack.enter_context(pysam.AlignmentFile(bam, 'rb')))

        if args.outfile and args.outfile != '-':
            logging.debug("Opening output file: {}".format(args.outfile))
            out_fh = stack.enter_context(open(args.outfile, 'w'))
        else:
            out_fh = sys.stdout

        writer = csv.DictWriter(out_fh, fieldnames=output_fields + samples)
        writer.writeheader()

        for record in SeqIO.parse(gb_fh, "genbank"):
            logging.info("Processing record ID: {}".format(record.id))

            for feature in record.features:

                if feature.type not in ['source', 'CDS', 'mat_peptide']:
                    continue

                [product] = feature.qualifiers.get('product', [None])
                if product is not None:
                    logging.debug("Fetching counts for product: {}".format(product))
                else:
                    logging.debug("Fetching counts for feature type: {}".format(feature.type))

                start = feature.location.start + 1
                end = feature.location.end

                region = "{}:{}-{}".format(record.id, start, end)
                logging.debug("Querying region: {}".format(region))

                counts = dict()

                for sample, bam in zip(samples, bams):
                    logging.debug("Getting counts for sample: {}".format(sample))

                    count = bam.count(region=region, read_callback=read_callback)
                    counts[sample] = count

                writer.writerow({
                    "region": region,
                    "strand": feature.strand,
                    "type": feature.type,
                    "product": product,
                    **counts})


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt as err:
        print("Interrupt signal received. Exiting...", file=sys.stderr)
        sys.exit(130)
