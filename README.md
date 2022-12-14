## Readme.md
## RepCounter, version 0.1

### General purpose of RepCounter
The purpose of RepCounter is to provide a lightweight, modular program for estimating the relative abundance of repeat families and classes in a next-generation sequencing project, such as RNA-seq or ChIP-seq. RepCounter is implemented in Python using Snakemake, and can be utilized on a cluster or on a local machine.
In it's intended implementation, RepCounter examines all annotated repeat types on your genome, rather than focusing on one or more select repeat elements. This allows you to have an unbiased view of the relative repeat expression in your NGS samples.

### Preparing annotation
The overall workflow of this program is to generate a custom genome from individual repeat elements on your genome of choice, such as mouse. Your sequencing reads for a given sample are then aligned to this repeat-based Repeatmasker genome. The sums of the individual repeat belonging to particular families or classes of repeats are then calculated for each sample, and normalized to the library size.
To begin, you will need to obtain the fasta file for your genome of choice, such as mm10, as well as the repeat masker annotation for the given genome. In the current implementation, we are parsing the UCSC repeatmasker track obtained from the UCSC table browser. Future implementations may allow for more flexible input. Please note that this program has been tested on mm10, and should be valid for most common genomes. However, should an issue arise with the parsing of the repeatmasker file, for example, due to changes in the header names, the Snakefile can be modified to suit your particular implementation.

To extract the necessary Repeatmasker file from UCSC, follow the following steps: UCSC table browser > choose organism > choose genome build > choose group "Variation and repeats" > choose "Repeatmasker" track > output "all fields from selected table"

## Dependency:
- [Anaconda](https://conda.io/docs/user-guide/install/linux.html) 

## Installation:
Clone this repository and change into the cloned NGS-Data-Charmer directory. 

To create an environment using the environment.yaml file, type the following:

`conda env create -f environment.yaml`

This will create a conda environment called ngs_data_charmer.

## Usage note:

You must manually activate the conda environment prior to running the sh files. Type the following to activate the environment:

`conda activate repcounter`

The reason for this requirement is a failure of the conda environment to successfully activate from within a shell script.

## Usage on a local machine:

Copy the config.yaml, run\_snakemake.sh and Snakefile to your NGS project directory. This directory should also contain a directory called 'fastq' wherein all the fastq files are placed. Make sure the project directory structure is as follows:
```
.
????????? config.yaml
????????? fastq
???   ????????? D1-WC_S2_L003_R1_001.fastq.gz
???   ????????? D1-WC_S2_L003_R2_001.fastq.gz
????????? run_snakemake.sh
????????? Snakefile
```
Make the required changes to the config.yaml file.

Finally, type `sh run_snakemake.sh` followed by the maximum number of CPU cores to be used by snakemake. For example, type `sh run_snakemake.sh 2` for 2 CPU cores. You can also type `nohup sh run_snakemake.sh 2 &` to run the pipeline in background.

## Usage on an LSF cluser:

Copy the config.yaml, run\_snakemake\_cluster.sh, cluster.json and Snakefile to your NGS project directory. This directory should also contain a directory called 'fastq' wherein all fastq files are placed. Make sure the project directory structure is as follows:
```
.
????????? cluster.json
????????? config.yaml
????????? fastq
???   ????????? negD1-WC-40_S2_L003_R1_001.fastq.gz
???   ????????? negD1-WC-40_S2_L003_R2_001.fastq.gz
????????? run_snakemake_cluster.sh
????????? Snakefile (required for all analyses)
```
Make the required changes to the config.yaml and cluster.json file. 

Optionally, to utilize multiple cores on the cluster enviroment for computationally heavy tasks (such as alignment), you may change the number of cores utilized in the cluster.json file. Following is an example requesting 4 compute cores on a single node, asking 12GB of memory per core -

```
    "__default__" :
    {  
        "queue"     : "premium",
        "allocation": "acc_labname",
        "n"         : 4,
        "resources" : "\"rusage[mem=12000] span[ptile=4]\"",
        "jobname"      : "{rule}.{wildcards}",
        "output"    : "logs/{rule}.{wildcards}.o",
        "error"     : "logs/{rule}.{wildcards}.e",
        "walltime"    : "02:00"
    }

```
You can make the above changes either to the '\__default\__' object alone or to any of the individual objects in the cluster.json file.

Finally, type `nohup sh run_snakemake_cluster.sh &` (to run in background).

### Program workflow
 ![ScreenShot](/dag/RNA.png)

### Output files
A compiled matrix of your sample's repeat expression, both raw counts and counts normalized to the sample size, are generated by RepCounter. Since each repeat belongs to a particular repeat name, repeat family, and repeat class, the expression of each of the three levels are saved for your usage.
#### Raw count files
1. counts_matrix.txt
2. Repeat_expression_repfamily_raw.txt
3. Repeat_expression_repname_raw.txt
4. Repeat_expression_repclass_raw.txt

#### Normalized count files
1. Repeat_expression_repfamily_norm.txt
1. Repeat_expression_repname_norm.txt
1. Repeat_expression_repclass_norm.txt

## Testing
For ease of testing, a subset of the RepeatMasker annotation for mouse (mm10) is provided with this repository (mm10_repeatmasker_ERVL.txt). The small size of this annotation should ensure a fast index build and sample alignment.
