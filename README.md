# Workflow for protein prediction #
This repository contains the analysis of protein prediction based on global similarity, domain, and structure.

## Requirments ##
*Python 3+
*cd-hit v4.8.1
*diamond
## Fetch protein sequences from NCBI using E-Utilities ##
```
usage: fetch_pro.py [-h] -i INPUT -o OUTPUT

A script to download protein sequences using NCBI E-Utilities.

options:
  -h, --help           show this help message and exit
  -i, --input INPUT    Input file containing protein IDs (one per line).
  -o, --output OUTPUT  Output file to save protein sequences.
```
## CD-hits to generate UniRef90 proteins ##
```
cd-hit -i input.faa -o output.faa -c 0.9 -aS 0.8 -aL 0.8

The protein DB was generated based on 90% identity and mutual 80% coverage (both query and subject)
```
