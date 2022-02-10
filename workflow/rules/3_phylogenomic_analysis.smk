
import os
import glob
import pandas as pd

###########################
rule dn_genomes:
    input:
        'status/gtdbtk.done',
        expand(os.path.join(DATA_DIR, "lists/{GENUS}/genomes_list.txt"), GENUS=GENUS_LIST)
    output:
        touch("status/dn_genomes.done")

########################################
# rules for ncbi-download-genomes #
########################################

rule download_genomes:
    input:
        expand(os.path.join(DATA_DIR, 'lists/{GENUS}/{GENUS}_accessions'),GENUS=GENUS_LIST)
    output:
        directory(expand(os.path.join(RESULTS_DIR,"Genomes/"),GENUS=GENUS_LIST))
    conda:
        os.path.join(ENV_DIR, "ncbi-g-d.yaml")
    script:
        os.path.join(SRC_DIR, "dn_genomes_script.py")