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

def config_with_default(key, default):
    try:
        return config[key]
    except KeyError:
        return default

# Location of temporary files
TMP_DIR = config["temp-dir"]

# Reports directory
REPORTS_DIR = TMP_DIR + "reports/"

# Input directory for fastq files
FASTQ_DIR = config["fastq-dir"]

# Trimmed reads directory
TRIMMED_READS_DIR = REPORTS_DIR + "trimmed_reads/"

# Max number of threads.
# TODO: is this needed?
MAX_THREADS = int(config_with_default("max-threads", "1"))

# Path to the trimmomatic jar
TRIM_PATH = config["trim-path"]

FASTQC_DIR = "{0}raw_reads/qc/".format(REPORTS_DIR)

#------------------------------------------------------------------------
# Build all the single filename wildcard patterns so we can use them multiple times
# without duplication.
# Note that we need to escape the Snakemake wildcards with double curly braces
# so that they pass through the string formatting.

# TODO: Can I ditch the use of separate sample and id wildcards? They only appear to be used together anyway.
FASTQ_FILE = "{0}{{sample}}_R{{id}}.fastq.gz".format(FASTQ_DIR)
RAW_FASTQC_ZIP = "{0}{{sample}}_R{{id}}_fastqc.zip".format(FASTQC_DIR)
RAW_FASTQC_HTML = "{0}{{sample}}_R{{id}}_fastqc.html".format(FASTQC_DIR)

# Patterns for forward and reverse fastq files used for the trim rule
FASTQ_FORWARD_FILE = "{0}{{sample}}_R1_fastqc.html".format(FASTQC_DIR)
FASTQ_REVERSE_FILE = "{0}{{sample}}_R2_fastqc.html".format(FASTQC_DIR)

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

#------------------------------------------------------------------------

# TODO - remove these
csiro_id = "ley015"
temp_loc = expand("/scratch1/{csiro_id}", csiro_id = csiro_id)
IDS = [1,2]
PUs = ['P','U']

sample = ["CA73YANXX_8_161220_BPO--000_Other_TAAGGCGA-CTCTCTAT_R_161128_SHADIL_LIB2500_M002"]#,
#            "CA73YANXX_8_161220_BPO--000_Other_TAAGGCGA-CTCTCTAT_R_161128_SHADIL_LIB2500_M002"]

#samples = {f[:-11] for f in os.listdir(".") if f.endswith("fastq.gz")}

trim_path = '/apps/trimmomatic/0.38/trimmomatic-0.38.jar'
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
        print(ALL_RAW_FASTQC_ZIP)
        print(ALL_RAW_FASTQC_HTML)

#------------------------------------------------------------------------

rule all:
    input:
        ALL_RAW_FASTQC_HTML,
        ALL_RAW_FASTQC_ZIP
        #expand("{temp_loc}/reports/raw_reads/qc/{sample}_R{id}_fastqc.zip", sample = sample, temp_loc = temp_loc, id = IDS),
        #expand("{temp_loc}/reports/raw_reads/qc/{sample}_R{id}_fastqc.html", sample = sample, temp_loc = temp_loc, id = IDS),
        #expand("{temp_loc}/reports/trimmed_reads/qc/{sample}_{id}{pu}_fastqc.zip", temp_loc = temp_loc, sample = sample, id = IDS, pu = PUs),
        #expand("{temp_loc}/reports/trimmed_reads/qc/{sample}_{id}{pu}_fastqc.html", temp_loc = temp_loc, sample = sample, id = IDS, pu = PUs),
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

rule trim:
    input:
        forward=FASTQ_FORWARD_FILE,
        reverse=FASTQ_REVERSE_FILE
    output:
        forward_paired="{0}{{sample}}_1P.fq.gz".format(TRIMMED_READS_DIR)
        forward_unpaired="{0}{{sample}}_1U.fq.gz".format(TRIMMED_READS_DIR)
        reverse_paired="{0}{{sample}}_2P.fq.gz".format(TRIMMED_READS_DIR)
        reverse_unpaired="{0}{{sample}}_2U.fq.gz".format(TRIMMED_READS_DIR)
        trimlog="{0}{{sample}}.log".format(TRIMMED_READS_DIR)
    threads:
        MAX_THREADS
    shell:
        """
        module load trimmomatic/0.38
        java -jar {TRIM_PATH} PE -phred33 \
        {input.forward} {input.reverse} \
        {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} \
        -trimlog {output.trimlog} \
        LEADING:3 \
        TRAILING:3 \
        ILLUMINACLIP:TrueSeq3-PE.fa:2:30:10
        module unload trimmomatic/0.38
        """

rule fastqc_trimmed:
    input:
        fq1 = expand("{temp_loc}/reports/trimmed_reads/{sample}_{id}{pu}.fq.gz", id = IDS, temp_loc = temp_loc, sample = sample, pu = PUs),
    output:
        zip2 = "{temp_loc}/reports/trimmed_reads/qc/{sample}_{id}{pu}_fastqc.zip",
        html2 = "{temp_loc}/reports/trimmed_reads/qc/{sample}_{id}{pu}_fastqc.html",
    threads:
        4
    shell:
        """
        module load fastqc/0.11.8
        mkdir -p {temp_loc}/reports/trimmed_reads/qc/
        fastqc -t 32 {input.fq1} -o {temp_loc}/reports/trimmed_reads/qc/
        module unload fastqc/0.11.8
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

        # """
        # /apps/trinity/2.3.2/util/align_and_estimate_abundance.pl --thread_count 16 --transcripts /OSM/CBR/AF/OZ_WHEAT/work/ref_seq/161010_Chinese_Spring_v1.0_pseudomolecules.fasta --seqType fq \
        # --left {input.left} --right {input.right} \
        # --est_method RSEM --output_dir {temp_loc}/scratch1/ley015/trimmed_reads --aln_method bowtie2 --max_ins_size 1500 \
        # --trinity_mode \
        # --prep_reference --output_prefix "RSEM_"
        # """

 #     --gene_trans_map /OSM/CBR/AF/OZ_WHEAT/work/ref_seq/161010_Chinese_Spring_v1.0_pseudomolecules_AGP.tsv \


        # """
        # module load trinity/2.8.4
        # Trinity --seqType fq --max_memory 50G \
        # --left {input.left} \
        # --right {input.right} \
        # --output {temp_loc}/reports/raw_reads/trinity/ \
        # ---CPU 6
        # module unload trinity/2.8.4
        # """

# rule trinity:
#     input:
#         left = expand("test_data/{sample}_R1.fastq.gz", sample = sample),
#         right = expand("test_data/{sample}_R1.fastq.gz", sample = sample),
#     output:
#         fa = expand("{temp_loc}/reports/raw_reads/trinity/Trinity.fasta", temp_loc = temp_loc),
#     shell:
#         """
#         module load trinity/2.8.4
#         Trinity --seqType fq --max_memory 50G \
#         --left {input.left} \
#         --right {input.right} \
#         --output {temp_loc}/reports/raw_reads/trinity/ \
#         ---CPU 6 --trimmomatic
#         module unload trinity/2.8.4
#         """
