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

OUTPUT_DIR = config["output-dir"]

# Reports directory
REPORTS_DIR = OUTPUT_DIR + "reports/"

FASTQC_DIR = REPORTS_DIR + "raw_reads/qc/"

# Trimmed reads directory
TRIMMED_READS_DIR = REPORTS_DIR + "trimmed_reads/"
TRIMMED_FASTQC_DIR = TRIMMED_READS_DIR + "qc/"

# Max number of threads.
# TODO: is this needed?
MAX_THREADS = int(config.get("max-threads", "1"))

# Trinity working directory
TRINITY_DIR = TRIMMED_READS_DIR + "RSEM_trinity/"

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

#-------------------
# trim rule patterns

FORWARD_PAIRED = "{0}{{sample}}_1P.fq.gz".format(TRIMMED_READS_DIR)
FORWARD_UNPAIRED = "{0}{{sample}}_1U.fq.gz".format(TRIMMED_READS_DIR)
REVERSE_PAIRED = "{0}{{sample}}_2P.fq.gz".format(TRIMMED_READS_DIR)
REVERSE_UNPAIRED = "{0}{{sample}}_2U.fq.gz".format(TRIMMED_READS_DIR)


#-------------------
# trinity_align rule patterns
FA_FILE = TRINITY_DIR + "{sample}/Trinity.fasta"

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

ALL_TRINITY_FILES = expand(FA_FILE, sample=samples)

#------------------------------------------------------------------------


rule test:
    run:
        from pprint import PrettyPrinter
        pp = PrettyPrinter(indent=1)
        print("OUTPUT_DIR: {0}".format(OUTPUT_DIR))
        print("MAX_THREADS: {0}".format(MAX_THREADS))
        print("FASTQ_DIR:", FASTQ_DIR)
        print("FASTQ_FILE: ", FASTQ_FILE)
        print("RAW_FASTQC_ZIP: ", RAW_FASTQC_ZIP)
        print("RAW_FASTQC_HTML: ", RAW_FASTQC_HTML)
        print()
        pp.pprint(ALL_RAW_FASTQC_FILES)
        print()
        pp.pprint(ALL_TRIMMED_FASTQC_FILES)
        print()
        pp.pprint(ALL_TRINITY_FILES)

#------------------------------------------------------------------------

rule all:
    input:
        ALL_RAW_FASTQC_FILES
        ,ALL_TRIMMED_FASTQC_FILES
        #,ALL_TRINITY_FILES

rule clean:
    shell: "rm -rf {OUTPUT_DIR}"

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
        # This is the FASTQ_FILE pattern, but forcing a matched pair
        "{0}{{sample}}_R1.fastq.gz".format(FASTQ_DIR),
        "{0}{{sample}}_R2.fastq.gz".format(FASTQ_DIR)
    output:
        forward_paired=FORWARD_PAIRED,
        forward_unpaired=FORWARD_UNPAIRED,
        reverse_paired=REVERSE_PAIRED,
        reverse_unpaired=REVERSE_UNPAIRED,
        trimlog="{0}{{sample}}.log".format(TRIMMED_READS_DIR)
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
        left=FORWARD_PAIRED,
        right=REVERSE_PAIRED
    output: FA_FILE
    threads: MAX_THREADS
    shell:
        """
        module load trinity/2.3.2

        # Map the insilico_read_normalization temp dir to a memory directory
        #mkdir ${TRINITY_DIR}
        #rm -f ${TRINITY_DIR}/read_partitions
        #mkdir $MEMDIR/read_partitions
        #ln -s $MEMDIR/read_partitions $OUTPUT_DIR

        # Run trinity
        Trinity --seqType fq --max_memory 50G --left {input.left} --right {input.right} --output {TRINITY_DIR} --CPU {threads}
        """
