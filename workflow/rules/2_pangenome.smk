import os
import glob
import pandas as pd


###########################
# default

rule pan_genome:
    input:
        'status/dn_genomes.done',
        expand(os.path.join(RESULTS_DIR, '{GENUS}/mmseqs2_all_seqs.fasta'), GENUS=GENUS_LIST)
    output:
        touch("status/pan_genome.done")

rule prokka_genomes:
    input:
        expand(os.path.join(RESULTS_DIR, '{GENUS}'), GENUS=GENUS_LIST)
    output:
        directory(os.path.join(RESULTS_DIR, '{GENUS}/Pangenome'))
    params:
        CPU=12
    conda:
        os.path.join(ENV_DIR, "prokka.yaml")
    script:
        os.path.join(SRC_DIR, "run_prokka.py")

rule mmseqs2:
    input:
        expand(os.path.join(RESULTS_DIR, '{GENUS}/Pangenome'), GENUS=GENUS_LIST)
    output:
        expand(os.path.join(RESULTS_DIR, '{GENUS}/mmseqs2_all_seqs.fasta'), GENUS=GENUS_LIST)
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    script:
        os.path.join(SRC_DIR, 'run_mmseqs2.py')

rule mOTUlizer:
    input:
        expand(os.path.join(RESULTS_DIR, '{GENUS}/Pangenome'), GENUS=GENUS_LIST)
    output:
        expand(os.path.join(RESULTS_DIR, '{GENUS}/mmseqs2_all_seqs.fasta'), GENUS=GENUS_LIST)
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    script:
        os.path.join(SRC_DIR, 'run_mOTUlizer.py')

rule mOTUpan:
    input:
        expand(os.path.join(RESULTS_DIR, '{GENUS}/Pangenome'), GENUS=GENUS_LIST)
    output:
        expand(os.path.join(RESULTS_DIR, '{GENUS}/mmseqs2_all_seqs.fasta'), GENUS=GENUS_LIST)
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    script:
        os.path.join(SRC_DIR, 'run_mOTUpanpy')
