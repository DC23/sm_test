# sm_test
Snakemake test files

Creating the snakefile to conduct the following steps:

2. Run FastQC
3. Trim with trimmomatic
4. align to reference and estimate abundance with Trinity - TBD
5 (?). Normalise FKPM and TPM - TBD

Additionals:
6. Map to ERCC reference in Bowtie2 - TBD
7. de-novo assembly with Trinity

# Execution Notes

* Clone from GitHub
* Copy `config-template.yaml` to `config.yaml` and modify to suit your data.
* Execute with `snakemake`. On Pearcey, snakemake is in the Python 3.6 and
  Miniconda3 modules.

