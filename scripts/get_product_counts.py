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

from collections import defaultdict
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
        '-o', '--outdir',
        help='Write output files to DIR',
        default=os.getcwd(),
        metavar='DIR')

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


def iter_qualifiers(gb_rec, name=None):
    for feature in gb_rec.features:
        qualifier_name = feature.qualifiers.get(name)
        if qualifier_name is not None:
            [qualifier] = qualifier_name
            yield qualifier, feature.location


def main():

    parser = get_option_parser()
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s: %(message)s', level=args.loglevel)

    os.makedirs(args.outdir, exist_ok=True)

    samples = args.bam_files

    # try to remove a common dirname
    if len(set(os.path.split(bam)[0] for bam in args.bam_files)) == 1:
        samples = [os.path.split(bam)[1] for bam in args.bam_files]
        # try to remove a common file extension
        if len(set(sample.split('.', 1)[1] for sample in samples)) == 1:
            samples = [sample.split('.', 1)[0] for sample in samples]

    count_types = [ "total_counts", "primary_counts", "unique_counts" ]
    output_fields = [ "record_id", "start", "end", "strand", "product" ]

    with ExitStack() as stack:

        logging.debug("Opening GenBank file: {}".format(args.gb_file))
        gb_fh = stack.enter_context(open(args.gb_file, 'rU'))

        bams = []

        for bam in args.bam_files:
            logging.debug("Opening BAM file: {}".format(bam))
            bams.append(stack.enter_context(pysam.AlignmentFile(bam, 'rb')))

        writers = []

        for count_type in count_types:
            gb_filestem = os.path.splitext(os.path.basename(args.gb_file))[0]
            count_filename = "{}.{}.csv".format(gb_filestem, count_type)
            count_filepath = os.path.join(args.outdir, count_filename)

            logging.debug("Opening counts file: {}".format(count_filepath))
            counts_fh = stack.enter_context(open(count_filepath, 'w'))

            writer = csv.DictWriter(counts_fh, fieldnames=output_fields + samples)
            writer.writeheader()

            writers.append(writer)

        for record in SeqIO.parse(gb_fh, "genbank"):
            logging.info("Processing record ID: {}".format(record.id))

            for product, location in iter_qualifiers(record, name='product'):
                logging.info("Fetching counts for product: {}".format(product))

                region = (record.id, location.start, location.end)
                logging.debug("Querying region: {}:{}-{}".format(*region))

                counts = defaultdict(dict)

                for sample, bam in zip(samples, bams):
                    logging.debug("Getting counts for sample: {}".format(sample))

                    total_count = 0
                    primary_count = 0
                    unique_count = 0

                    for alignment in bam.fetch(*region):
                        if not alignment.is_secondary:
                            primary_count += 1
                        if alignment.mapping_quality == 255:
                            unique_count += 1
                        total_count += 1

                    counts["total_counts"][sample] = total_count
                    counts["primary_counts"][sample] = primary_count
                    counts["unique_counts"][sample] = unique_count

                row_data = {
                    "record_id": record.id,
                    "start": location.start,
                    "end": location.end,
                    "strand": location.strand,
                    "product": product,
                }

                for count_type, writer in zip(count_types, writers):
                    writer.writerow({**row_data, **counts[count_type]})


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt as err:
        print("Interrupt signal received. Exiting...", file=sys.stderr)
        sys.exit(130)
