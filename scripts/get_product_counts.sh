#!/bin/bash
for i in M35027.1 KU321639.1 AM258992.1 ; do
  python scripts/get_product_counts.py -o "${i}.total_counts.csv" "${i}.gb" "$@"
  python scripts/get_product_counts.py -o "${i}.total_proper_pairs.csv" -f 1 -f 2 "${i}.gb" "$@"
  python scripts/get_product_counts.py -o "${i}.primary_counts.csv" -F 256 "${i}.gb" "$@"
  python scripts/get_product_counts.py -o "${i}.primary_proper_pairs.csv" -F 256 -f 1 -f 2 "${i}.gb" "$@"
  python scripts/get_product_counts.py -o "${i}.unique_counts.csv" -q 255 "${i}.gb" "$@"
  python scripts/get_product_counts.py -o "${i}.unique_proper_pairs.csv" -q 255 -f 1 -f 2 "${i}.gb" "$@"
done
