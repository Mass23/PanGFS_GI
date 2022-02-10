import os
import glob
import pandas as pd

###########################
rule clean_mags:
    input:
        os.path.join(RESULTS_DIR, "checkm_out")
    output:
        touch("status/1_clean_mags.done")

################################################
# rule for collecting mags with taxo. matching #
################################################

rule run_gtdbtk:
    input:
        os.path.join(MAGS_DIR)
    output:
        directory(os.path.join(RESULTS_DIR, 'gtdbtk_output'))
    log:
        os.path.join(RESULTS_DIR, 'logs/gtdbtk.log')
    conda:
        os.path.join(ENV_DIR, 'gtdbtk_updated.yaml')
    params:
        config["gtdbtk"]["path"]
    threads:
        config["gtdbtk"]["threads"]
    message:
        "Running GTDB toolkit on MAGs"
    shell:
        "export OMP_NUM_THREADS=48 && "
        "export PENBLAS_NUM_THREADS=48 && "
        "export MKL_NUM_THREADS=48 && "
        "export VECLIB_MAXIMUM_THREADS=48 && "
        "export NUMEXPR_NUM_THREADS=48 && "
        "(date && export GTDBTK_DATA_PATH={params} && gtdbtk classify_wf --cpus {threads} -x fasta --genome_dir {input} --out_dir {output} && date) &> {log}"

rule list_target_mags:
    input:
        GTDBTK_DIR=os.path.join(RESULTS_DIR, 'gtdbtk_output/')
    params:
        GENUS=GENUS_LIST
    output:
        expand(os.path.join(DATA_DIR, 'accessions/{GENUS}/mags_list.txt'),GENUS=GENUS_LIST)
    run:
        for i in range(0,len(GENUS_LIST)):
            tax_string = 'g__' + params.GENUS[i]
            gtdbtk_file = pd.read_csv(f'{input.GTDBTK_DIR}/gtdbtk.bac120.summary.tsv', sep='\t')
            gtdbtk_sub = gtdbtk_file[gtdbtk_file['classification'].str.contains(tax_string)].user_genome
            gtdbtk_sub.to_csv(output[i], header =  False, sep='\t', index=False)

checkpoint copy_target_mags:
    input:
        expand(os.path.join(DATA_DIR, 'accessions/{GENUS}/mags_list.txt'),GENUS=GENUS_LIST)
    output:
        os.path.join(RESULTS_DIR, 'MAGs')
    run:
        os.mkdir(output)
        for i in range(0,len(input)):
            files_to_move = [i for i in open(input[i], 'r').read().split('\n') if i != '']
            for file_to_move in files_to_move:
                contigs_length = [len(rec.seq) for rec in SeqIO.parse(os.path.join(MAGS_DIR, file_to_move) + '.fasta','fasta')]
                if sum(contigs_length) > 100000:
                    os.system('cp -v ' + os.path.join(MAGS_DIR, file_to_move) + '.fasta ' + output])

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.copy_target_mags.get(**wildcards).output[0]
    return [file for file in glob.glob(os.path.join(RESULTS_DIR, 'MAGs/*.fasta')) if os.path.isfile(file)]

rule mag_purify:
    input:
        aggregate_input
    output:
        directory(os.path.join(RESULTS_DIR, "cleaned_MAGs"))
    wildcard_constraints:
        files="|".join([file for file in glob.glob(os.path.join(RESULTS_DIR, 'MAGs/*.fasta'))])
    conda:
        os.path.join(ENV_DIR, "magpurify.yaml")
    shell:
        "magpurify phylo-markers --db /mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0 {input} $(basename {input})_magpurify && "
        "magpurify clade-markers --db /mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0 {input} $(basename {input})_magpurify && "
        "magpurify tetra-freq {input} $(basename {input})_magpurify && "
        "magpurify gc-content {input} $(basename {input})_magpurify && "
        "magpurify known-contam --db /mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0 {input} $(basename {input})_magpurify && "
        "magpurify clean-bin {input} $(basename {input})_magpurify $(echo {input} | sed -e 's|/MAGS/|/cleaned_MAGs/|g' -e 's|.fasta|_clean.fasta|g') "

rule checkm:
    input:
        os.path.join(RESULTS_DIR, "cleaned_MAGs/")
    output:
        directory(os.path.join(RESULTS_DIR, "checkm_out"))
    threads:
        config["checkm"]["threads"]
    conda:
        os.path.join(ENV_DIR, "checkm.yaml")
    shell:
        "checkm lineage_wf -r -t threads -x fa {input} {output}"

