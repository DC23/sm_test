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

FASTQC_DIR = REPORTS_DIR + "raw_reads/qc/"

# Trimmed reads directory
TRIMMED_READS_DIR = REPORTS_DIR + "trimmed_reads/"
TRIMMED_FASTQC_DIR = TRIMMED_READS_DIR + "qc/"

# Max number of threads.
# TODO: is this needed?
MAX_THREADS = int(config.get("max-threads", "1"))

# Path to the trimmomatic jar
TRIM_PATH = config["trim-path"]

#------------------------------------------------------------------------
# Build the single filename wildcard patterns so we can use them multiple times
# without duplication.

# If you use string concatenation then wildcards should use single curly braces
FASTQ_FILE = FASTQ_DIR + "{sample}_R{id}.fastq.gz"

# If using string formatting, then double curly braces are needed on wildcards
RAW_FASTQC_ZIP = "{0}{{sample}}_R{{id}}_fastqc.zip".format(FASTQC_DIR)
RAW_FASTQC_HTML = "{0}{{sample}}_R{{id}}_fastqc.html".format(FASTQC_DIR)

#-------------------
# trimqc rule output patterns
# Defined here so they can be used in the rule and to build all-file lists
TRIMMED_FASTQC_ZIP = "{0}{{sample}}_{{id}}{{pu}}_fastqc.zip".format(TRIMMED_FASTQC_DIR)
TRIMMED_FASTQC_HTML = "{0}{{sample}}_{{id}}{{pu}}_fastqc.html".format(TRIMMED_FASTQC_DIR)

#------------------------------------------------------------------------
# Build the lists of all expected output files.
# These are written to globals to avoid duplication when the lists are used in
# multiple rules.
# Note the common pattern of building the all-file lists from the single file patterns

samples, ids = glob_wildcards(FASTQ_FILE)
# convert to sets to remove duplication
samples = set(samples)
ids = set(ids)

ALL_RAW_FASTQC_FILES = \
    expand(RAW_FASTQC_ZIP, sample=samples, id=ids) + \
    expand(RAW_FASTQC_HTML, sample=samples, id=ids)

ALL_TRIMMED_FASTQC_FILES = \
    expand(TRIMMED_FASTQC_ZIP, sample=samples, id=(1,2), pu=("P", "U")) + \
    expand(TRIMMED_FASTQC_HTML, sample=samples, id=(1,2), pu=("P", "U"))

#------------------------------------------------------------------------

# TODO - remove these once the trinity rule is rewritten
csiro_id = "ley015"
temp_loc = expand("/scratch1/{csiro_id}", csiro_id = csiro_id)
IDS = [1,2]
PUs = ['P','U']
sample = ["CA73YANXX_8_161220_BPO--000_Other_TAAGGCGA-CTCTCTAT_R_161128_SHADIL_LIB2500_M002"]#,

#------------------------------------------------------------------------

rule test:
    run:
        from pprint import PrettyPrinter
        pp = PrettyPrinter(indent=1)
        print("TMP_DIR: {0}".format(TMP_DIR))
        print("MAX_THREADS: {0}".format(MAX_THREADS))
        print("FASTQ_DIR:", FASTQ_DIR)
        print("FASTQ_FILE: ", FASTQ_FILE)
        print("RAW_FASTQC_ZIP: ", RAW_FASTQC_ZIP)
        print("RAW_FASTQC_HTML: ", RAW_FASTQC_HTML)
        print()
        pp.pprint(ALL_RAW_FASTQC_FILES)
        print()
        pp.pprint(ALL_TRIMMED_FASTQC_FILES)

#------------------------------------------------------------------------

rule all:
    input:
        ALL_RAW_FASTQC_FILES
        ,ALL_TRIMMED_FASTQC_FILES
        #expand("{temp_loc}/reports/trimmed_reads/RSEM_trinity/{sample}/Trinity.fasta", temp_loc = temp_loc, sample = sample),

rule clean:
    shell: "rm -rf {TMP_DIR}"

# A single input file, with wildcards for sample and ID
rule fastqc_raw:
    input: FASTQ_FILE
    output: RAW_FASTQC_ZIP, RAW_FASTQC_HTML
    threads: 1  # we are only processing one file per rule execution
    shell:
        """
        module load fastqc/0.11.8
        fastqc -t {threads} {input} -o {FASTQC_DIR}
        """

rule trim:
    input:
        "{0}{{sample}}_R1.fastq.gz".format(FASTQ_DIR),
        "{0}{{sample}}_R2.fastq.gz".format(FASTQ_DIR)
    output:
        files=expand("{path}{{sample}}_{id}{pu}.fq.gz", path=TRIMMED_READS_DIR, id=(1,2), pu=("P", "U")),
        trimlog="{0}{{sample}}.log".format(TRIMMED_READS_DIR)
    threads: 2  # 2 input files per call
    shell:
        """
        module load trimmomatic/0.38
        trimmomatic PE -phred33 \
        -threads {threads} \
        {input} \
        {output.files} \
        -trimlog {output.trimlog} \
        LEADING:3 \
        TRAILING:3 \
        ILLUMINACLIP:TrueSeq3-PE.fa:2:30:10
        """

rule fastqc_trimmed:
    input: "{0}{{sample}}_{{id}}{{pu}}.fq.gz".format(TRIMMED_READS_DIR)
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
