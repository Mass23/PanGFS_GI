import os
import glob
import pandas as pd

###########################
rule clean_mags:
    input:
        os.path.join(DATA_DIR, 'status/cleaning.done')
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
            tax_string = f'g__{params.GENUS[i]}'
            gtdbtk_file = pd.read_csv(f'{input.GTDBTK_DIR}/gtdbtk.bac120.summary.tsv', sep='\t')
            gtdbtk_sub = gtdbtk_file[gtdbtk_file['classification'].str.contains(tax_string)].user_genome
            gtdbtk_sub.to_csv(output[i], header =  False, sep='\t', index=False)

checkpoint copy_target_mags:
    input:
        expand(os.path.join(DATA_DIR, 'accessions/{GENUS}/mags_list.txt'),GENUS=GENUS_LIST)
    output:
        target=directory(os.path.join(RESULTS_DIR, "MAGs"))
    run:
        os.mkdir(os.path.join(RESULTS_DIR,'MAGs'))
        for i in range(0,len(input)):
            files_to_move = [i for i in open(input[i], 'r').read().split('\n') if i != '']
            for file_to_move in files_to_move:
                contigs_length = [len(rec.seq) for rec in SeqIO.parse(f'{os.path.join(MAGS_DIR, file_to_move)}.fasta','fasta')]
                if sum(contigs_length) > 100000:
                    os.system(f'cp -v {os.path.join(MAGS_DIR, file_to_move)}.fasta {output}')

rule bin_link:
    input:
        os.path.join(RESULTS_DIR, "MAGs/{i}.fasta")
    output:
        os.path.join(RESULTS_DIR, "linked/MAGs/{i}.fasta")
    wildcard_constraints:
        GENUS="|".join(GENUS_LIST)
    shell:
        "ln -vs {input} {output}"


rule mag_purify:
    input:
        os.path.join(RESULTS_DIR, "linked/MAGs/{i}.fasta")
    output:
        os.path.join(RESULTS_DIR, "cleaned_MAGs/{i}_clean.fasta")
    conda:
        os.path.join(ENV_DIR, "magpurify.yaml")
    shell:
        "magpurify phylo-markers --db /mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0 {input} $(basename {input})_magpurify && "
        "magpurify clade-markers --db /mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0 {input} $(basename {input})_magpurify && "
        "magpurify tetra-freq {input} $(basename {input})_magpurify && "
        "magpurify gc-content {input} $(basename {input})_magpurify && "
        "magpurify known-contam --db /mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0 {input} $(basename {input})_magpurify && "
        "magpurify clean-bin {input} {input}_magpurify {output}"

def aggregate_mags(wildcards):
    checkpoint_output = checkpoints.copy_target_mags.get(**wildcards).output.target
    return expand(os.path.join(RESULTS_DIR, "cleaned_MAGs/{i}_clean.fasta"), i=glob_wildcards(os.path.join(RESULTS_DIR,"MAGs/{i}.fasta")).i)


rule cleaning_done:
    input:
        aggregate_mags
    output:
        touch(os.path.join(DATA_DIR, 'status/cleaning.done'))
