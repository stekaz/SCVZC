#!/bin/bash
set -euo pipefail

bam_file="${1}"
regions_file="${2}"

readarray -t regions < "${regions_file}"

# total counts
echo -n "total_counts:"
samtools view -c "${bam_file}"
echo -n "total_proper_pairs:"
samtools view -c -f 1 -f 2 "${bam_file}"
echo -n "total_region_counts:"
samtools view -c "${bam_file}" "${regions[@]}"
echo -n "total_region_proper_pairs:"
samtools view -c -f 1 -f 2 "${bam_file}" "${regions[@]}"

# primary counts
echo -n "primary_counts:"
samtools view -F 256 -c "${bam_file}"
echo -n "primary_proper_pairs:"
samtools view -F 256 -c -f 1 -f 2 "${bam_file}"
echo -n "primary_region_counts:"
samtools view -F 256 -c "${bam_file}" "${regions[@]}"
echo -n "primary_region_proper_pairs:"
samtools view -F 256 -c -f 1 -f 2 "${bam_file}" "${regions[@]}"

# unique counts
echo -n "unique_counts:"
samtools view -q 255 -c "${bam_file}"
echo -n "unique_proper_pairs:"
samtools view -q 255 -c -f 1 -f 2 "${bam_file}"
echo -n "unique_region_counts:"
samtools view -q 255 -c "${bam_file}" "${regions[@]}"
echo -n "unique_region_proper_pairs:"
samtools view -q 255 -c -f 1 -f 2 "${bam_file}" "${regions[@]}"
