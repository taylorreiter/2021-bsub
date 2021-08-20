import pandas as pd
import csv

m = pd.read_csv("inputs/working_metadata.tsv", sep = "\t", header = 0)
SAMPLES = m.sort_values(by='read_count')['run_accession']
LIBRARIES = m['library_name'].unique().tolist()
STUDY = m['study_accession'].unique().tolist()

class Checkpoint_GatherResults:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint (e.g. solve for {acc} wildcard). This approach
    is documented at this url: 
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_genome_accs(self):
        gather_csv = f'outputs/genbank/bsub_assemblies.x.genbank.gather.csv'
        assert os.path.exists(gather_csv)

        genome_accs = []
        with open(gather_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               acc = row['name'].split(' ')[0]
               genome_accs.append(acc)
        print(f'loaded {len(genome_accs)} accessions from {gather_csv}.')

        return genome_accs

    def __call__(self, w):
        global checkpoints

        # wait for the results of rule 'format_bsub_accessions'; 
        # this will trigger exception until that rule has been run.
        checkpoints.format_bsub_accessions.get(**w)

        # parse accessions in gather output file
        genome_accs = self.get_genome_accs()

        p = expand(self.pattern, acc=genome_accs, **w)
        return p


rule all:
    input:
         Checkpoint_GatherResults("outputs/charcoal/{acc}_genomic.fna.gz.clean.fa.gz")

########################################
## PREPROCESSING
########################################

rule download_fastq_files_R1:
    output: 
        r1="inputs/raw/{sample}_1.fastq.gz",
    threads: 1
    resources:
        mem_mb=1000
    run:
        row = m.loc[m['run_accession'] == wildcards.sample]
        fastq_1 = row['fastq_ftp_1'].values
        fastq_1 = fastq_1[0]
        shell("wget -O {output.r1} {fastq_1}")


rule download_fastq_files_R2:
    output:
        r2="inputs/raw/{sample}_2.fastq.gz"
    threads: 1
    resources:
        mem_mb=1000
    run:
        row = m.loc[m['run_accession'] == wildcards.sample]
        fastq_2 = row['fastq_ftp_2'].values
        fastq_2 = fastq_2[0]
        shell("wget -O {output.r2} {fastq_2}")


rule cat_libraries_R1:
    input: expand("inputs/raw/{sample}_1.fastq.gz", sample = SAMPLES)
    output: expand("inputs/cat/{library}_1.fastq.gz", library = LIBRARIES)
    threads: 1
    resources:
        mem_mb=4000
    run: 
        merge_df = m[['library_name','run_accession']]
        merge_df = copy.deepcopy(merge_df)
        merge_df['run_accession'] = merge_df['run_accession'].apply(lambda x: f"inputs/raw/{x}_1.fastq.gz")
        merge_dict = merge_df.groupby('library_name')['run_accession'].apply(lambda g: g.values.tolist()).to_dict()
        for library in merge_dict.keys():
            # merge SRR files
            to_merge = merge_dict[library]
            # Check if the merged file results from a single or multiple fastq files.
            # For n-to-1 merging, concatenate input files to produce the output file
            merge_nb = len(to_merge)
            if merge_nb > 1:
                cmd = "cat " + " ".join(to_merge) + " > " + "inputs/cat/" + library + "_1.fastq.gz"
            else:
                cmd = "ln --relative --force -s " + " ".join(to_merge) + " inputs/cat/" + library + "_1.fastq.gz"
            os.system(cmd)
    
rule cat_libraries_R2:
    input: expand("inputs/raw/{sample}_2.fastq.gz", sample = SAMPLES)
    output: expand("inputs/cat/{library}_2.fastq.gz", library = LIBRARIES)
    threads: 1
    resources:
        mem_mb=4000
    run: 
        merge_df = m[['library_name','run_accession']]
        merge_df = copy.deepcopy(merge_df)
        merge_df['run_accession'] = merge_df['run_accession'].apply(lambda x: f"inputs/raw/{x}_2.fastq.gz")
        merge_dict = merge_df.groupby('library_name')['run_accession'].apply(lambda g: g.values.tolist()).to_dict()
        for library in merge_dict.keys():
            # merge SRR files
            to_merge = merge_dict[library]
            # Check if the merged file results from a single or multiple fastq files.
            # For n-to-1 merging, concatenate input files to produce the output file
            merge_nb = len(to_merge)
            if merge_nb > 1:
                cmd = "cat " + " ".join(to_merge) + " > " + "inputs/cat/" + library + "_2.fastq.gz"
            else:
                cmd = "ln --relative --force -s " + " ".join(to_merge) + " inputs/cat/" + library + "_2.fastq.gz"
            os.system(cmd)
    
rule adapter_trim_files:
    input:
        r1 = "inputs/cat/{library}_1.fastq.gz",
        r2 = 'inputs/cat/{library}_2.fastq.gz',
        adapters = 'inputs/adapters2.fa'
    output:
        r1 = 'outputs/trim/{library}_R1.trim.fq.gz',
        r2 = 'outputs/trim/{library}_R2.trim.fq.gz',
        o1 = 'outputs/trim/{library}_o1.trim.fq.gz',
        o2 = 'outputs/trim/{library}_o2.trim.fq.gz'
    conda: 'envs/trim.yml'
    threads: 1
    resources:
        mem_mb=8000
    shell:'''
     trimmomatic PE {input.r1} {input.r2} \
             {output.r1} {output.o1} {output.r2} {output.o2} \
             ILLUMINACLIP:{input.adapters}:2:0:15 MINLEN:31  \
             LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2
    '''

rule cutadapt_files:
    input:
        r1 = 'outputs/trim/{library}_R1.trim.fq.gz',
        r2 = 'outputs/trim/{library}_R2.trim.fq.gz',
    output:
        r1 = 'outputs/cut/{library}_R1.cut.fq.gz',
        r2 = 'outputs/cut/{library}_R2.cut.fq.gz',
    conda: 'envs/trim2.yml'
    threads: 1
    resources:
        mem_mb=8000
    shell:'''
    cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o {output.r1} -p {output.r2} {input.r1} {input.r2}
    '''

rule fastqc:
    input:
        r1 = 'outputs/cut/{library}_R1.cut.fq.gz',
        r2 = 'outputs/cut/{library}_R2.cut.fq.gz',
    output:
        r1 = 'outputs/fastqc/{library}_R1.cut_fastqc.html',
        r2 = 'outputs/fastqc/{library}_R2.cut_fastqc.html'
    conda: 'envs/trim2.yml'
    threads: 1
    resources:
        mem_mb=4000
    shell:'''
    fastqc -o outputs/fastqc {input} 
    '''
    
rule remove_host:
# http://seqanswers.com/forums/archive/index.php/t-42552.html
# https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?usp=sharing
    output:
        r1 = 'outputs/bbduk/{library}_R1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{library}_R2.nohost.fq.gz',
        human_r1='outputs/bbduk/{library}_R1.human.fq.gz',
        human_r2='outputs/bbduk/{library}_R2.human.fq.gz'
    input: 
        r1 = 'outputs/cut/{library}_R1.cut.fq.gz',
        r2 = 'outputs/cut/{library}_R2.cut.fq.gz',
        human='inputs/host/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz'
    threads: 1
    resources:
        mem_mb=64000
    conda: 'envs/trim.yml'
    shell:'''
    bbduk.sh -Xmx64g t=3 in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.human_r1} outm2={output.human_r2} k=31 ref={input.human}
    '''

rule kmer_trim_reads:
    input: 
        'outputs/bbduk/{library}_R1.nohost.fq.gz',
        'outputs/bbduk/{library}_R2.nohost.fq.gz'
    output: "outputs/abundtrim/{library}.abundtrim.fq.gz"
    conda: 'envs/trim.yml'
    threads: 1
    resources:
        mem_mb=64000
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 60e9 -V - -o {output}
    '''

#############################################
# B. subtilis reference preprocessing
#############################################

rule grab_bsub_accessions:
    input: "/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.v2.csv"
    output: "outputs/genbank/bsub_acc.csv",
    resources:
        mem_mb = 4000
    threads: 1
    shell:"""
    grep -i s__Bacillus {input} | grep -i subtilis > {output}
    """

checkpoint format_bsub_accessions:
    input: "outputs/genbank/bsub_acc.csv",
    output: "outputs/genbank/bsub_assemblies.x.genbank.gather.csv",
    resources:
        mem_mb = 4000
    threads: 1
    shell:"""
    sed '1iname,lineage' {input} > {output}
    """

rule generate_bsub_lineages:
    input: "outputs/genbank/bsub_acc.csv"
    output: "outputs/genbank/bsub_assemblies.x.genbank.lineages.csv",
    resources:
        mem_mb = 4000
    threads: 1
    shell:"""
    sed 's/,/_genomic.fna.gz,/1' {input} > {output}
    """
    
# specifying make_sgc_conf as the target downloads the genomes of interest,
# circumventing a bug in genome grist that prevents using the target
# download gather genomes. The sgc conf file is a dummy file -- it will be
# written to outputs/sgc, but the conf file has the wrong catlas bases.
rule download_bsub_assemblies:
    input: 
        gather_grist = "outputs/genbank/bsub_assemblies.x.genbank.gather.csv",
        conf = "conf//genome-grist-conf.yml"
    output: "genbank_genomes/{acc}_genomic.fna.gz"
    conda: "envs/genome-grist.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell:'''
    genome-grist run {input.conf} --until make_sgc_conf --nolock
    touch {output}
    '''

rule generate_charcoal_genome_list:
    input:  ancient(Checkpoint_GatherResults("genbank_genomes/{acc}_genomic.fna.gz"))
    output: "outputs/charcoal_conf/charcoal.genome-list.txt"
    threads: 1
    resources:
        mem_mb=500
    shell:'''
    ls genbank_genomes/*gz | xargs -n 1 basename > {output} 
    '''

rule charcoal_decontaminate_bsub:
    input:
        genomes = ancient(Checkpoint_GatherResults("genbank_genomes/{acc}_genomic.fna.gz")),
        genome_list = "outputs/charcoal_conf/charcoal.genome-list.txt",
        conf = "conf/charcoal-conf.yml",
        genome_lineages = "outputs/genbank/bsub_assemblies.x.genbank.lineages.csv",
        db="/group/ctbrowngrp/gtdb/databases/gtdb-rs202.genomic.k31.zip",
        db_lineages="/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.csv"
    output: 
        hitlist="outputs/charcoal/stage1_hitlist.csv",
        clean_finished="outputs/charcoal/clean_finished.txt"
    resources:
        mem_mb = 128000
    threads: 8
    conda: "envs/charcoal.yml"
    shell:'''
    python -m charcoal run {input.conf} -j {threads} clean --nolock --latency-wait 15 --rerun-incomplete
    touch {output.clean_finished}
    '''

rule touch_decontaminated_bsub:
    input: 
        "outputs/charcoal/stage1_hitlist.csv",
        "outputs/charcoal/clean_finished.txt"
    output: "outputs/charcoal/{acc}_genomic.fna.gz.clean.fa.gz"
    shell:'''
    touch {output}
    '''
#############################################
# Spacegraphcats Genome Queries
#############################################

    
rule make_sgc_conf_files:
    input:
        csv = "outputs/genbank/bsub_assemblies.x.genbank.gather.csv",
        queries = ancient(Checkpoint_GatherResults("outputs/charcoal/{acc}_genomic.fna.gz.clean.fa.gz")),
    output:
        conf = "outputs/sgc_conf/{library}_k31_r1_conf.yml"
    resources:
        mem_mb = 500
    threads: 1
    run:
        query_list = "\n- ".join(input.queries)
        with open(output.conf, 'wt') as fp:
           print(f"""\
catlas_base: {wildcards.library}
input_sequences:
- outputs/abundtrim/{wildcards.library}.abundtrim.fq.gz
ksize: 31
radius: 1
paired_reads: true
search:
- {query_list}
""", file=fp)

rule spacegraphcats_genome_nbhd_query:
    input: 
        queries = ancient(Checkpoint_GatherResults("outputs/charcoal/{acc}_genomic.fna.gz.clean.fa.gz")), 
        conf = ancient("outputs/sgc_conf/{library}_k31_r1_conf.yml"),
        reads = "outputs/abundtrim/{library}.abundtrim.fq.gz"
    output:
        "outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/results.csv"
    params: outdir = "outputs/sgc_genome_queries"
    conda: "envs/spacegraphcats.yml"
    resources:
        mem_mb = 150000
    threads: 1
    shell:'''
    python -m spacegraphcats run {input.conf} extract_contigs extract_reads --nolock --outdir={params.outdir}  
    '''

#rule prokka_queries:
#    output: 
#        ffn = 'outputs/prokka_shared_assemblies/{acc}.ffn',
#        faa = 'outputs/prokka_shared_assemblies/{acc}.faa'
#    input: 'outputs/charcoal/{acc}_genomic.fna.gz.clean.fa.gz'
#    conda: 'envs/prokka.yml'
#    resources:
#        mem_mb = 8000
#    threads: 2
#    params: 
#        outdir = 'outputs/prokka_shared_assemblies/',
#        prefix = lambda wildcards: wildcards.acc,
#        gzip = lambda wildcards: "outputs/charcoal/" + wildcards.acc + "_genomic.fna.gz.clean.fa.gz"
#    shell:'''
#    gunzip {input}
#    prokka {params.gzip} --outdir {params.outdir} --prefix {params.prefix} --metagenome --force --locustag {params.prefix} --cpus {threads} --centre X --compliant
#    mv {params.prefix}.ffn {output.ffn}
#    mv {params.prefix}.faa {output.faa}
#    gzip {params.gzip}
#    '''
