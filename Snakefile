# To run:
# module load miniconda3/4.3.24
# snakemake -j 999 --cluster-config cluster-configs/pearcey.json --cluster 'sbatch --job-name {cluster.job-name} --ntasks-per-node {cluster.ntasks-per-node} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}'

# To dryrun:
# module load miniconda3/4.3.24
# snakemake -n

# To create a DAG file in dryrun:
# module load miniconda3/4.3.24
# snakemake -n --dag | dot -Tsvg > dag.svg


#------------------------------------------------------------------------
# Load the configuration
configfile: "config.yaml"

# Input directory for fastq files
FASTQ_DIR = config["fastq-dir"]

# Location of temporary files
# TODO: It looks more like the main output directory
TMP_DIR = config["temp-dir"]

# Reports directory
REPORTS_DIR = TMP_DIR + "reports/"

FASTQC_DIR = "{0}raw_reads/qc/".format(REPORTS_DIR)

# Trimmed reads directory
TRIMMED_READS_DIR = REPORTS_DIR + "trimmed_reads/"
TRIMMED_FASTQC_DIR = TRIMMED_READS_DIR + "qc/"

# Max number of threads.
# TODO: is this needed?
MAX_THREADS = int(config.get("max-threads", "1"))

# Path to the trimmomatic jar
TRIM_PATH = config["trim-path"]

#------------------------------------------------------------------------
# Build all the single filename wildcard patterns so we can use them multiple times
# without duplication.
# Note that we need to escape the Snakemake wildcards with double curly braces
# so that they pass through the string formatting.

FASTQ_FILE = "{0}{{sample}}_R{{id}}.fastq.gz".format(FASTQ_DIR)
RAW_FASTQC_ZIP = "{0}{{sample}}_R{{id}}_fastqc.zip".format(FASTQC_DIR)
RAW_FASTQC_HTML = "{0}{{sample}}_R{{id}}_fastqc.html".format(FASTQC_DIR)

#-------------------
# trim rule
# input patterns
FASTQ_FORWARD_FILE = "{0}{{sample}}_R1.fastq.gz".format(FASTQ_DIR)
FASTQ_REVERSE_FILE = "{0}{{sample}}_R2.fastq.gz".format(FASTQ_DIR)

# output patterns
FORWARD_PAIRED = "{0}{{sample}}_1P.fq.gz".format(TRIMMED_READS_DIR)
FORWARD_UNPAIRED = "{0}{{sample}}_1U.fq.gz".format(TRIMMED_READS_DIR)
REVERSE_PAIRED = "{0}{{sample}}_2P.fq.gz".format(TRIMMED_READS_DIR)
REVERSE_UNPAIRED = "{0}{{sample}}_2U.fq.gz".format(TRIMMED_READS_DIR)
TRIMLOG = "{0}{{sample}}.log".format(TRIMMED_READS_DIR)

#-------------------
# trimqc rule
# matches all outputs from trim rule
TRIMMED_FILE = "{0}{{sample}}_{{id}}{{pu}}.fq.gz".format(TRIMMED_READS_DIR)

# output patterns
TRIMMED_FASTQC_ZIP = "{0}{{sample}}_{{id}}{{pu}}_fastqc.zip".format(TRIMMED_FASTQC_DIR)
TRIMMED_FASTQC_HTML = "{0}{{sample}}_{{id}}{{pu}}_fastqc.html".format(TRIMMED_FASTQC_DIR)

#------------------------------------------------------------------------
# Build the lists of all expected output files.
# These are written to globals to avoid duplication when the lists are used in
# multiple rules.

samples, ids = glob_wildcards(FASTQ_FILE)
# convert to sets to remove duplication
samples = set(samples)
ids = set(ids)

ALL_RAW_FASTQC_ZIP = expand(RAW_FASTQC_ZIP, sample=samples, id=ids)
ALL_RAW_FASTQC_HTML = expand(RAW_FASTQC_HTML, sample=samples, id=ids)
ALL_RAW_FASTQC_FILES = ALL_RAW_FASTQC_ZIP + ALL_RAW_FASTQC_HTML

# Don't need this, but it shows the common pattern of building the all file
# lists from the single file patterns
#ALL_TRIMMED_FILES = \
    #expand(FORWARD_PAIRED, sample=samples) + \
    #expand(FORWARD_UNPAIRED, sample=samples) + \
    #expand(REVERSE_PAIRED, sample=samples) + \
    #expand(REVERSE_UNPAIRED, sample=samples) + \
    #expand(TRIMLOG, sample=samples)

ALL_TRIMMED_FASTQC_FILES = \
    expand(TRIMMED_FASTQC_ZIP, sample=samples, id=(1,2), pu=("P", "U")) + \
    expand(TRIMMED_FASTQC_HTML, sample=samples, id=(1,2), pu=("P", "U"))

#------------------------------------------------------------------------

# TODO - remove these
csiro_id = "ley015"
temp_loc = expand("/scratch1/{csiro_id}", csiro_id = csiro_id)
IDS = [1,2]
PUs = ['P','U']

sample = ["CA73YANXX_8_161220_BPO--000_Other_TAAGGCGA-CTCTCTAT_R_161128_SHADIL_LIB2500_M002"]#,


#------------------------------------------------------------------------

rule test:
    run:
        print("TMP_DIR: {0}".format(TMP_DIR))
        print("MAX_THREADS: {0}".format(MAX_THREADS))
        print("FASTQ_DIR:", FASTQ_DIR)
        print("FASTQ_FILE: ", FASTQ_FILE)
        print("RAW_FASTQC_ZIP: ", RAW_FASTQC_ZIP)
        print("RAW_FASTQC_HTML: ", RAW_FASTQC_HTML)
        print(samples)
        print(ids)
        print(ALL_RAW_FASTQC_FILES)
        print()
        print(ALL_TRIMMED_FASTQC_FILES)

#------------------------------------------------------------------------

rule all:
    input:
        ALL_RAW_FASTQC_FILES
        #,ALL_TRIMMED_FILES
        ,ALL_TRIMMED_FASTQC_FILES
        #expand("{temp_loc}/reports/trimmed_reads/RSEM_trinity/{sample}/Trinity.fasta", temp_loc = temp_loc, sample = sample),

rule clean:
    shell:
        """
        rm -rf {TMP_DIR}
        rm -f {ALL_RAW_FASTQC_HTML}
        rm -f {ALL_RAW_FASTQC_ZIP}
        """

# A single input file, with wildcards for sample and ID
rule fastqc_raw:
    input: FASTQ_FILE
    output:
        zip=RAW_FASTQC_ZIP,
        html=RAW_FASTQC_HTML
    threads: 1  # we are only processing one file per rule execution
    shell:
        """
        module load fastqc/0.11.8
        fastqc -t {threads} {input} -o {FASTQC_DIR}
        """

# TODO: can this rule be broken down further so that a single call
# handles just one input file (either forward or reverse) to produce a single
# output pair of files (paired, unpaired)?
rule trim:
    input:
        forward=FASTQ_FORWARD_FILE,
        reverse=FASTQ_REVERSE_FILE
    output:
        forward_paired=FORWARD_PAIRED,
        forward_unpaired=FORWARD_UNPAIRED,
        reverse_paired=REVERSE_PAIRED,
        reverse_unpaired=REVERSE_UNPAIRED,
        trimlog=TRIMLOG
    threads: 2  # 2 input files per call
    shell:
        """
        module load trimmomatic/0.38
        trimmomatic PE -phred33 \
        -threads {threads} \
        {input} \
        {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} \
        -trimlog {output.trimlog} \
        LEADING:3 \
        TRAILING:3 \
        ILLUMINACLIP:TrueSeq3-PE.fa:2:30:10
        """

rule fastqc_trimmed:
    input: TRIMMED_FILE
    output: TRIMMED_FASTQC_ZIP, TRIMMED_FASTQC_HTML
    threads: 1  # one file per execution
    shell:
        """
        module load fastqc/0.11.8
        fastqc -t {threads} {input} -o {TRIMMED_FASTQC_DIR}
        """

rule trinity_align:
    input:
        left = expand("{temp_loc}/reports/trimmed_reads/{sample}_1P.fq.gz", temp_loc = temp_loc, sample = sample),
        right = expand("{temp_loc}/reports/trimmed_reads/{sample}_2P.fq.gz", temp_loc = temp_loc, sample = sample),
    output:
        fa = expand("{temp_loc}/reports/trimmed_reads/RSEM_trinity/{sample}/Trinity.fasta", temp_loc = temp_loc, sample = sample),
    shell:
        """
        module load trinity/2.3.2
        Trinity --seqType fq --max_memory 50G --left {input.left} --right {input.right} --output /scratch1/ley015/reports/trimmed_reads/RSEM_trinity --CPU 6
        module unload trinity/2.3.2
        """
