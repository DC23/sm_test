# To run:
# module load miniconda3/4.3.24
# snakemake -j 999 --cluster-config cluster-configs/pearcey.json --cluster 'sbatch --job-name {cluster.job-name} --ntasks-per-node {cluster.ntasks-per-node} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}'

# To dryrun:
# module load miniconda3/4.3.24
# snakemake -n

# To create a DAG file in dryrun:
# module load miniconda3/4.3.24
# snakemake -n --dag | dot -Tsvg > dag.svg 


csiro_id = "ley015"
temp_loc = expand("/scratch1/{csiro_id}", csiro_id = csiro_id)
IDS = [1,2]
PUs = ['P','U']

sample = ["CA73YANXX_8_161220_BPO--000_Other_TAAGGCGA-CTCTCTAT_R_161128_SHADIL_LIB2500_M002"]#,
#            "CA73YANXX_8_161220_BPO--000_Other_TAAGGCGA-CTCTCTAT_R_161128_SHADIL_LIB2500_M002"]

#samples = {f[:-11] for f in os.listdir(".") if f.endswith("fastq.gz")}

MAX_THREADS = 32
trim_path = '/apps/trimmomatic/0.38/trimmomatic-0.38.jar'

rule all:
    input:
        expand("{temp_loc}/reports/raw_reads/qc/{sample}_R{id}_fastqc.zip", sample = sample, temp_loc = temp_loc, id = IDS),
        expand("{temp_loc}/reports/raw_reads/qc/{sample}_R{id}_fastqc.html", sample = sample, temp_loc = temp_loc, id = IDS),
        expand("{temp_loc}/reports/trimmed_reads/qc/{sample}_{id}{pu}_fastqc.zip", temp_loc = temp_loc, sample = sample, id = IDS, pu = PUs),
        expand("{temp_loc}/reports/trimmed_reads/qc/{sample}_{id}{pu}_fastqc.html", temp_loc = temp_loc, sample = sample, id = IDS, pu = PUs),
        expand("{temp_loc}/reports/trimmed_reads/RSEM_trinity/{sample}/Trinity.fasta", temp_loc = temp_loc, sample = sample),

# rule clean:
#    shell:
#       """
#       rm -rf blahblahblah
#       """     


rule fastqc_raw:
    input:
        fastq = expand("test_data/{sample}_R{id}.fastq.gz", sample = sample, id = IDS),
    output:
        zip1 = expand("{temp_loc}/reports/raw_reads/qc/{sample}_R{id}_fastqc.zip", sample = sample, temp_loc = temp_loc, id = IDS),
        html1 = expand("{temp_loc}/reports/raw_reads/qc/{sample}_R{id}_fastqc.html", sample = sample, temp_loc = temp_loc, id = IDS),
    threads:
        MAX_THREADS
    shell:
        """
        module load fastqc/0.11.8
        mkdir -p {temp_loc}/reports/raw_reads/
        fastqc -t {threads} {input.fastq} -o {temp_loc}/reports/raw_reads/qc/
        """
 
rule trim:
    input:
        forward = "test_data/{sample}_R1.fastq.gz",
        reverse = "test_data/{sample}_R2.fastq.gz",
    output:
        forward_paired = "{temp_loc}/reports/trimmed_reads/{sample}_1P.fq.gz",
        forward_unpaired = "{temp_loc}/reports/trimmed_reads/{sample}_1U.fq.gz",
        reverse_paired = "{temp_loc}/reports/trimmed_reads/{sample}_2P.fq.gz",
        reverse_unpaired = "{temp_loc}/reports/trimmed_reads/{sample}_2U.fq.gz",
        trimlog = "{temp_loc}/reports/trimmed_reads/{sample}.log",        
    threads:
        MAX_THREADS
    shell:
        """
        module load trimmomatic/0.38               
        java -jar {trim_path} PE -phred33 \
        {input.forward} {input.reverse} \
        {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} \
        -trimlog {temp_loc}/reports/trimmed_reads/{sample}.log \
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