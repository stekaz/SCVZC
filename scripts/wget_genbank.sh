#!/bin/bash
ncbi_viewer="https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi"
for acc in "$@"; do
  wget -O "${acc}.gb" "${ncbi_viewer}?db=nuccore&id=${acc}&retmode=txt"
done
