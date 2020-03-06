#!/bin/bash
set -e

gencode_mouse="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24"
wget "${gencode_mouse}/GRCm38.primary_assembly.genome.fa.gz"
wget "${gencode_mouse}/gencode.vM24.primary_assembly.annotation.gtf.gz"

zcat "GRCm38.primary_assembly.genome.fa.gz" > GRCm38_SCVZC.fa

# VACV, ZIKV, CHIKV
for acc in "M35027" "KU321639" "AM258992" ; do

    # wget the viral FASTA files
    ncbi_viewer="https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi"
    wget -O "${acc}.fa" "${ncbi_viewer}?db=nuccore&id=${acc}&rettype=fasta"

    # append formatted header
    sed -n '1s/^>\([^ ]*\).*/>\1 \1/p' "${acc}.fa" >> GRCm38_SCVZC.fa

    # append folded sequence
    sed '1d' "${acc}.fa" | tr -d '\n' | fold -w 60 >> GRCm38_SCVZC.fa

    # append newline character
    echo >> GRCm38_SCVZC.fa

    # clean up    
    rm -f "${acc}.fa"
done

gzip GRCm38_SCVZC.fa
