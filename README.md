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

## Visualising the DAG

On CSIRO clusters, you need to load the imagemagick module to get the `display`
command. Once done, this will work:

`snakemake --dag | dot -Tpng | display`

## Optional: Create Python 3.7 Virtual Environment on Pearcey

This is a workaround for the fact that the Python 3.7 module does not contain
snakemake, and the version installed in Python 3.6 is old.

Load the Python 3.7 module and create a Python virtual environment:

* `module load git python/3.7.2`
* `git clone git@github.com:EPLeyne/sm_test.git`
* `cd sm_test`
* `source /apps/python/3.7.2/bin/virtualenvwrapper_lazy.sh`
* `mkvirtualenv -a . snakemake`
* `workon snakemake`
* `pip install snakemake`

From that point on you will need to activate the virtual Python environment
each time:

* `module load git python/3.7.2`
* `source /apps/python/3.7.2/bin/virtualenvwrapper_lazy.sh`
* `workon snakemake`
