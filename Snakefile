"""
Author: Y. Ahmed-Braimah

"""

import json
import os
import re
import pandas as pd
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output

##--------------------------------------------------------------------------------------##
## Functions
##--------------------------------------------------------------------------------------##

# To print process messages
def message(x):
  print()

# To remove suffix from a string
def rstrip(text, suffix):
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

## define environment variables

##--------------------------------------------------------------------------------------##
## Global config files:
##--------------------------------------------------------------------------------------##

configfile: 'config.yml'

QUERY_GTF = config['QUERY_GTF']
QUERY_PEP = config['QUERY_PEP']
REF_FOLD = config['REF_FOLD']
ANNOTS = config['ANNOTS']
GENE_LIST = config['GENE_LIST']

units_table = pd.read_table(GENE_LIST)
SAMPLES = list(units_table.iloc[:, 0].unique())

# Full path to a folder where final output files will be deposited.
OUT_DIR = config['OUT_DIR']

## Create the final output directory if it doesn't already exist
if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Final expected output(s)
rule all:
    input:
        expand(join(OUT_DIR, '{sample}', '{sample}.done'), sample = SAMPLES)

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule BlastSyn:
    params:
        gtf = QUERY_GTF,
        pep = QUERY_PEP,
        ref = REF_FOLD,
        ant = ANNOTS
    output:
        done = join(OUT_DIR, '{sample}', '{sample}.done')
    log:
        join(OUT_DIR, '{sample}', '{sample}.stdout.log')
    benchmark:
        join(OUT_DIR, '{sample}', '{sample}.benchmark.tsv')
    message:
        """--- Finding syntenic orthologs for gene "{wildcards.sample}" """
    run:
        shell('BlastSynteny'
                ' -g {params.gtf}'
                ' -p {params.pep}'
                ' -f {params.ref}'
                ' -a {params.ant}'
                ' -i {wildcards.sample}'
                ' > {log} 2>&1')
        shell('mv .tmp_synblast_folder.{wildcards.sample} {wildcards.sample}.synteny_blast_results')
        shell('mv {wildcards.sample}.synteny_blast_results ' + join(OUT_DIR, '{wildcards.sample}'))
        shell('touch {output.done}')


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##
