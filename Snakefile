## Snakefile
configfile: "config.yaml"
myfastqpath = "fastq/"

##Import necessary modules
import os
import sys
from os import listdir
from os.path import isfile, join
import re
from Bio import SeqIO
import gzip
from collections import Counter


## Define Helper functions
# Create the pattern of file endings
def create_endings(x):
    """
    Returns a list of likely fastq file endings

    Input Parameter:
    x (int): Either 1 or 2, indicating the forward (1) or reverse (2) read.

    Returns:
    list: A list of strings, representing the file endings a user might
    use for denoting their fastq files.
    """
    return(["_R" + str(x) + "_001.fastq", "_R" + str(x) + "_001.fq",
            "_R" + str(x) + ".fastq", "_R" + str(x) + ".fq",
            "_" + str(x) + ".fastq", "_" + str(x) + ".fq",
            ".R" + str(x) + "_001.fastq", ".R" + str(x) + "_001.fq",
            ".R" + str(x) + ".fastq", ".R" + str(x) + ".fq",
            "." + str(x) + ".fastq", "." + str(x) + ".fq",
            "_r" + str(x) + "_001.fastq", "_r" + str(x) + "_001.fq",
            "_r" + str(x) + ".fastq", "_r" + str(x) + ".fq",
            ".r" + str(x) + "_001.fastq", ".r" + str(x) + "_001.fq",
            ".r" + str(x) + ".fastq", ".r" + str(x) + ".fq"])

# Function to list the fastq files present in the fastq folder
def getfilelist(myfastqpath):
    """
    Extracts fastq files from the files present in your fastq directory.

    Input Parameter:
    myfastqpath (string): directory containing your fastq files.

    Returns:
    list: List containing two strings.
    1st string is all non-metadata files in the fastq directory
    2nd string is all non-metadata files ending in '.gz'
    """
    onlyfiles = [f for f in listdir(myfastqpath) if
                 isfile(join(myfastqpath, f))]
    onlyfiles = [i for i in onlyfiles if
                 i.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))]
    gzfiles = [i for i in onlyfiles if i.endswith((".gz"))]
    return([onlyfiles, gzfiles])

## Function to rename files
def rename_files(oldname, replacement, myfastqpath):
    [os.rename(join(myfastqpath, i), join(myfastqpath, y))
     for i, y in zip(oldname, replacement)]

# Function to unify fastq files to single file ending
def fix_input_files(file_suffix, input_fileset, myfastqpath):
    """
    Renames mixed input fastq files to the most common file ending and
    returns the selected file ending. NOTE: This step permenantly
    renames your fastq files from their original file ending.

    Input Parameter:
    file_suffix (string): ".gz" or ""; Gzipped fastq files are expected
    to end with the suffix ".gz". If files are NOT gzipped, the input is "".

    input_fileset (list): List of fastq file names to be examined.
    As written, gzipped files are listed within the variable 'gzfiles'
    and non-gzipped files are listed within the variable 'onlyfiles'.

    Returns:
    list: A list containing four strings, the selected Read1 (forward read)
    file ending and the corresponding Read2 (reverse read) file ending,
    a list of all fastq-like files, and a list of gzipped fastq-like files.
    """

    # file_suffix, input_fileset = [".gz", gzfiles]
    # Create the series of fastq file endings to search
    base_endings_r1, base_endings_r2 = [create_endings(i) for i in (1, 2)]
    # Define the R1 and R2 suffix pairs for reference
    ending_dictionary = dict(zip(base_endings_r1, base_endings_r2))
    mylist = list()  # Create empty list

    # Traverse the R1 base endings to find the common ending
    for x in base_endings_r1:
        matched_ends = [
            i for i in input_fileset if i.endswith(x + file_suffix)]
        if(len(matched_ends) > 0):
            mylist.extend([x]*len(matched_ends))

    # If all samples are single-end
    if len(mylist) == 0:
        print("Your dataset appears to be entirely single-end files.")
        odd_files = [i for i in input_fileset
                     if i.endswith(".fq" + file_suffix)]
        if len(odd_files) > 0:
            old_rep = [i.replace(".fq" + suffix, ".fastq" + suffix)
                       for i in odd_files]
            rename_files(odd_files, old_rep, myfastqpath)

        # Re-assess fastq directory content and return filenames
        return([".fastq", ".fastq", getfilelist(myfastqpath)[0],
            getfilelist(myfastqpath)[1]])
    # If R1 endings are present, check values and correct file names
    else:
        # create dictionary of mixed file endings
        mylist_endings = list(Counter(mylist).keys())
        # Find most common file ending
        myR1_suffix = max(Counter(mylist).items(), key=lambda x: x[1])[0]
        # Match chosen R1 ending to correct R2 ending
        myR2_suffix = ending_dictionary[myR1_suffix]
        # remove main R1 suffix from dictionary
        mylist_endings.remove(myR1_suffix)

        # Process forward reads
        if len(mylist_endings) > 0:
            for x in mylist_endings:
                oldnames = [
                    i for i in input_fileset if i.endswith(x + file_suffix)]
                old_rep = [i.replace(x, myR1_suffix) for i in oldnames]
                rename_files(oldnames, old_rep, myfastqpath)

            mylist = list()  # Create empty list to hold R2 file endings
            # Traverse the R2 base endings to endings
            for x in base_endings_r2:
                matched_ends = [i for i in input_fileset if i.endswith(
                    x + file_suffix) and x != myR2_suffix]
                if(len(matched_ends) > 0):
                    mylist.append(x)  # Create list of R2 files to be renamed
            if len(mylist) > 0: # Rename R2 files that don't match desired ending
                for x in mylist:
                    oldnames = [
                        i for i in input_fileset if i.endswith(x + file_suffix)]
                    old_rep = [i.replace(x, myR2_suffix) for i in oldnames]
                    rename_files(oldnames, old_rep, myfastqpath)

        # Re-assess file names
        if file_suffix == ".gz":
            input_fileset = getfilelist(myfastqpath)[1]
        else:
            input_fileset = getfilelist(myfastqpath)[0]

        # Now process single end files
        # Identify files that do not match the current R1, R2 ending
        odd_files = [i for i in input_fileset if not
                     i.endswith(myR1_suffix + file_suffix) if not
                     i.endswith(myR2_suffix + file_suffix)]

        # Partition single end files according to ending
        fastq_odd_1 = [i for i in odd_files if i.endswith(
            ".fastq" + file_suffix)]
        fastq_odd_2 = [i for i in odd_files if i.endswith(".fq" + file_suffix)]

        # If any apparently single-end files exist, then rename them
        if len(odd_files) > 0:
            print("Now unifying " + str(len(odd_files)) +
                  " single-end files to \"" + myR1_suffix +
                  file_suffix + "\" ending")
            # rename 'fastq' single-end files to correct ending
            if len(fastq_odd_1) > 0:
                old_rep = [i.replace(".fastq" + file_suffix,
                            myR1_suffix + file_suffix) for i in fastq_odd_1]
                rename_files(fastq_odd_1, old_rep, myfastqpath)
            # rename 'fq' single-end files to correct ending
            if len(fastq_odd_2) > 0:
                old_rep = [i.replace(".fq" + file_suffix,
                            myR1_suffix + file_suffix) for i in fastq_odd_2]
                rename_files(fastq_odd_2, old_rep, myfastqpath)
        # Re-assess and return filenames and file endings
        return([myR1_suffix, myR2_suffix, getfilelist(myfastqpath)[0],
            getfilelist(myfastqpath)[1]])

# Create read-pair inputs for sample processing
def create_fastq_inputs(config):
    """
    Creates the fastq file inputs needed for read trimming steps of
    the snakemake pipeline

    Input Parameter:
    config (dict): Dictionary derived from config.yaml and any
    additional key:value pairs added during the file preperation steps.

    Returns:
    list: List of two strings;
    1st string denotes the forward read
    2nd string denotes the reverse read
    """
    return([os.path.join("fastq","{sample}"+expand("{ending}{suffix}", \
        ending=R1_file_ending, suffix=suffix)[0]+""),
        os.path.join("fastq","{sample}"+expand("{ending}{suffix}", \
            ending=R2_file_ending, suffix=suffix)[0]+"")])

#### End helper functions

# Retrieve the list of fastq files
onlyfiles, gzfiles = getfilelist(myfastqpath)

# Raise exception if no fastq files present
if len(gzfiles) == 0 and len(onlyfiles) == 0:
    raise NameError(
        "You do not seem to have any fastq files present to process. Exiting.")

# Raise exception if fastq files are a mixture of gzipped and non-gzipped files
if len(gzfiles) > 0 and len(gzfiles) != len(onlyfiles):
    myinput = "You have a mixture of gzipped files and non-gzipped files\n \
                Only {} of total {} files are gzipped!"
    raise NameError(print(myinput.format(len(gzfiles), len(onlyfiles))))

# Unify fastq file endings and return the final ending to be used.
if len(gzfiles) > 0:
    R1_file_ending, R2_file_ending, onlyfiles, gzfiles = \
            fix_input_files(".gz", gzfiles, myfastqpath)
    suffix = ".gz"
else:
    R1_file_ending, R2_file_ending, onlyfiles, gzfiles = \
            fix_input_files("", onlyfiles, myfastqpath)
    suffix = ""

sample_string = os.path.join("fastq","{sample}"+R1_file_ending+suffix)
SAMPLES, = glob_wildcards(sample_string)

# Check the file pairing
# Raise exception for non-paired PE files
if config["type"] == "single":
    print("You selected single-end reads\nRead pairing not being checked...")
elif config["type"] == "paired":
    len_r1 = len([i for i in onlyfiles if i.endswith(R1_file_ending + suffix)])
    if len_r1*2 != len(onlyfiles):
        myinput = "One or more samples do not have a read pair!\nIf using \
            paired-end samples, please ensure each sample has read 1 and \
            read 2 files\nAborting..."
        raise NameError(myinput)  # Raise exception to break workflow
else:
    myinput = "You have specified unknown read type: " + \
        config["type"] + "\nPlease specify either \"paired\" or \"single\" \
        in the config.yaml file, then rerun the pipeline."
    raise NameError(myinput)

##### RULES SECTION
# Generate input rule for Snakemake
rule all:
    input:
        ["output/logs/genome_finalcheck.txt",
        "output/repeatmasker.gtf",
        "output/Repeat_expression_repclass_norm.txt"]

## Step 1: Prepare the bed file from UCSC repeatmasker
rule prep_bed:
    input:
        config["annotation_file"]
    output:
        bed="output/repeatmasker.bed",
        gtf="output/repeatmasker.gtf"
    log:
        "output/logs/prep_bed.log"
    run:
        shell("awk '{{print $6,$7,$8,$11,$2,$10}}' {input} | sed 's/ /\\t/g' | sed 1d > {output.bed}")
        shell("awk 'BEGIN{{OFS=\"\\t\"}}{{print $11\"::\"$6\":\"$7\"-\"$8,\"repeatmasker\",\"gene\", 1, $8-$7,\".\",\"+\",\".\",\"gene_id \"$11\"::\"$6\":\"$7\"-\"$8}}' {input} | sed 1d > {output.gtf}")

rule generate_fasta:
    input:
        "output/repeatmasker.bed"
    params:
        fasta = config["genome_fasta"],
        index_folder = config["index_folder"]
    output:
        "output/repeatmasker.fa"
    run:
        shell("bedtools getfasta -name -fi {params.fasta} -bed {input} > {output}")

rule build_index:
    input:
        ancient("output/repeatmasker.fa")
    output:
        "output/logs/genome_finalcheck.txt"
    params:
        index_folder = config["index_folder"]
    threads: 
        config["threads"]
    log:
        "output/logs/build_index.log"
    run:
        shell("mkdir -p {params.index_folder}")
        shell("hisat2-build -p {threads} {input} {params.index_folder}/genome")
        shell("touch {output}")

## Select correct rule(s) for trimming reads
rule trim_fastq_fastqc:
    input:
        pair1 = create_fastq_inputs(config)[0],
        pair2 = create_fastq_inputs(config)[1]
    output:
        trimmed_pair1 = temp("output/trim_fastq/{sample}_R1_trimmed.fq.gz"),
        trimmed_pair2 = temp("output/trim_fastq/{sample}_R2_trimmed.fq.gz"),
        fastqc_zipfile1 = "output/fastqc/{sample}_R1_fastqc.zip",
        fastqc_zipfile2 = "output/fastqc/{sample}_R2_fastqc.zip"
    log:
        "output/logs/{sample}.trim_adapters.log"
    params:
        pair2 = create_fastq_inputs(config)[1]
    run:
        shell("mkdir -p output/temp_dir")
        if config["type"] == "paired":
            # mv files to R1 and R2 ending in temporary directory
            shell("cp {input.pair1} \
                output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
            shell("cp {params.pair2} \
                output/temp_dir/{wildcards.sample}_R2.fq{suffix}")
            shell("trim_galore \
                --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                output/temp_dir/{wildcards.sample}_R2.fq{suffix} --paired --trim-n \
                -o ./output/trim_fastq")
            shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                output/temp_dir/{wildcards.sample}_R2.fq{suffix} \
                -o ./output/fastqc")
            shell("mv output/trim_fastq/{wildcards.sample}_R1_val_1.fq.gz \
                output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz"),
            shell("mv output/trim_fastq/{wildcards.sample}_R2_val_2.fq.gz \
                output/trim_fastq/{wildcards.sample}_R2_trimmed.fq.gz")
        if config["type"] == "single":
            # mv files to R1 and R2 ending in temporary directory
            shell("cp {input.pair1} \
                output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
            shell("trim_galore --trim-n \
                --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                -o ./output/trim_fastq --basename {wildcards.sample}")
            shell("mv output/trim_fastq/{wildcards.sample}_trimmed.fq.gz \
                output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz")
            shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                -o ./output/fastqc")
            shell("touch {output.trimmed_pair2}")
            shell("touch {output.fastqc_zipfile2}")


# hisat2index = config["index_folder"],
# max_k = config["max_k"]

rule fastq_to_bam_HISAT:
    input:
        trimmed_pair = ["output/trim_fastq/{sample}_R1_trimmed.fq.gz", \
                        "output/trim_fastq/{sample}_R2_trimmed.fq.gz"],
        genome_done = "output/logs/genome_finalcheck.txt"
    params:
        hisat2index = config["index_folder"],
        max_k = config["max_k"]
    output:
        r1 = "output/bam/{sample}_R1_mapped.bam",
        r2 = "output/bam/{sample}_R2_mapped.bam",
        r1_full = "output/bam/{sample}_R1_sorted.bam"
    threads: config["threads"]
    log:
        "output/logs/{sample}.alignment.log"
    run:
        ## Align forward read 
        shell("hisat2 -k {params.max_k} -p {threads} --repeat --no-spliced-alignment \
            -x {params.hisat2index}/genome \
            -U {input.trimmed_pair[0]} > output/bam/{wildcards.sample}_R1.sam 2> {log}")
        ## Store alignments as sorted bam files for long-term storage
        shell("samtools view -b output/bam/{wildcards.sample}_R1.sam > output/bam/{wildcards.sample}_R1.bam")
        shell("samtools sort -@ 8 -O BAM -o output/bam/{wildcards.sample}_R1_sorted.bam output/bam/{wildcards.sample}_R1.bam")
        shell("samtools index output/bam/{wildcards.sample}_R1_sorted.bam")
        ## Subset the sam files to mapped reads only
        shell("samtools view -O BAM -F 4 output/bam/{wildcards.sample}_R1_sorted.bam > output/bam/{wildcards.sample}_R1_mapped.bam")
        shell("samtools index output/bam/{wildcards.sample}_R1_mapped.bam")
        ## Clean up unneeded files
        shell("rm output/bam/{wildcards.sample}_R1.sam output/bam/{wildcards.sample}_R1.bam")
        if config["type"] == "single":
            shell("touch {output.r2}")        
        if config["type"] == "paired":
            ## Align reverse read 
            shell("hisat2 -k {params.max_k} -p {threads} --repeat --no-spliced-alignment \
                -x {params.hisat2index}/genome \
                -U {input.trimmed_pair[1]} > output/bam/{wildcards.sample}_R2.sam 2> {log}")
            ## Store alignments as sorted bam files for long-term storage
            shell("samtools view -b output/bam/{wildcards.sample}_R2.sam > output/bam/{wildcards.sample}_R2.bam")
            shell("samtools sort -@ 8 -O BAM -o output/bam/{wildcards.sample}_R2_sorted.bam output/bam/{wildcards.sample}_R2.bam")
            shell("samtools index output/bam/{wildcards.sample}_R2_sorted.bam")
            ## Subset the sam files to mapped reads only
            shell("samtools view -O BAM -F 4 output/bam/{wildcards.sample}_R2_sorted.bam > output/bam/{wildcards.sample}_R2_mapped.bam")
            shell("samtools index output/bam/{wildcards.sample}_R2_mapped.bam")
            ## Remove original sam files
            shell("rm output/bam/{wildcards.sample}_R2.sam output/bam/{wildcards.sample}_R2.bam")

rule count_reads:
    input:
        r1 = "output/bam/{sample}_R1_mapped.bam",
        r2 = "output/bam/{sample}_R2_mapped.bam"
    output:
        r1 = "output/counts/{sample}_R1_htseq.tsv",
        r2 = "output/counts/{sample}_R2_htseq.tsv"
    params:
        gtf="output/repeatmasker.gtf"
    log:
        "output/logs/{sample}_htseq.log"
    run:
        shell("featureCounts -M -O --fraction  -t gene -a {params.gtf} \
            -o {output.r1} {input.r1} 2> output/logs/{wildcards.sample}_R1.feature_counts.log")
        if config["type"] == "paired":
            shell("featureCounts -M -O --fraction  -t gene -a {params.gtf} \
            -o {output.r2} {input.r2} 2> output/logs/{wildcards.sample}_R2.feature_counts.log")
        else:
            shell("touch {output.r2}")

rule counts_matrix:
    input:
        hts_ct = expand("output/counts/{sample}_{readtype}_htseq.tsv", sample=SAMPLES,readtype=["R1","R2"]) if config["type"] == "paired" else expand("output/counts/{sample}_R1_htseq.tsv", sample=SAMPLES)
    output:
        total_matrix="output/counts_matrix.txt"
    run:
        import pandas as pd
        import platform
        ## To extract sample name, recurse over R1 list
        R1_list=[x for x in input.hts_ct if "_R1_htseq.tsv" in x]
        R2_list=[x for x in input.hts_ct if "_R2_htseq.tsv" in x]
        sample_base=[x.replace("_R1_htseq.tsv","") for x in [x.replace("output/counts/","") for x in R1_list]]
        ##Compile R1 matrix
        dict_of_counts = {}
        for file in R1_list:
            sample = file.replace("_R1_htseq.tsv","")
            if platform.system() != 'Windows':
                sample = sample.split("/")[2]
            else:
                sample = sample.split("\\")[2]
            dict_of_counts[sample] = {}
            with open(file, "r") as infile:
                next(infile)
                next(infile)
                for lines in infile:
                    lines = lines.strip().split("\t")
                    if not lines[0].startswith("__"):
                        dict_of_counts[sample][lines[0]] = float(lines[-1])
        ## Transform R1 counts into a dataframe
        r1_df = pd.DataFrame(dict_of_counts)
        print(r1_df)
        if config["type"] == "paired": 
            dict_of_counts = {}
            for file in R2_list:
                sample = file.replace("_R2_htseq.tsv","")
                if platform.system() != 'Windows':
                    sample = sample.split("/")[2]
                else:
                    sample = sample.split("\\")[2]
                dict_of_counts[sample] = {}
                with open(file, "r") as infile:
                    next(infile)
                    next(infile)
                    for lines in infile:
                        lines = lines.strip().split("\t")
                        if not lines[0].startswith("__"):
                            dict_of_counts[sample][lines[0]] = float(lines[-1])
            r2_df = pd.DataFrame(dict_of_counts)
            print(r2_df)
            tot_df=r1_df+r2_df
            tot_df.to_csv(output[0], sep='\t')
        else:
            r1_df.to_csv(output[0], sep='\t')

rule normalize_counts:
    input:
        anno=config["annotation_file"],
        bamfiles=expand("output/logs/{sample}.alignment.log", sample=SAMPLES),
        total_matrix="output/counts_matrix.txt"
    output:
        family_file="output/Repeat_expression_repfamily_norm.txt",
        class_file="output/Repeat_expression_repclass_norm.txt",
        repname_file="output/Repeat_expression_repname_norm.txt"
    run:
        import pandas as pd
        ## Import the unmodified repeatmasker file
        orig=pd.read_csv(input.anno,sep='\t')
        orig=orig[['genoName','genoStart','genoEnd','repName','swScore','strand','repClass','repFamily']]
        orig['id']=orig['repName']+"::"+orig['genoName']+":"+orig['genoStart'].astype(str)+"-"+orig['genoEnd'].astype(str) 
        ## Import the counts matrix
        ct=pd.read_csv(input.total_matrix,index_col=0,sep='\t')
        ct.index.rename("id",inplace=True)
        sample_list=list(ct.columns)
        ##summarize the counts, according to repgroup
        ct_anno=ct.merge(orig,on='id')
        repfamily_raw=ct_anno.groupby("repFamily")[sample_list].sum()
        repname_raw=ct_anno.groupby("repName")[sample_list].sum()
        repclass_raw=ct_anno.groupby("repClass")[sample_list].sum()
        ## Read in the R1 samfile to get the total 
        bamsize=dict()
        for i in sample_list:
            mybam=pd.read_csv("output/logs/"+i+".alignment.log",header=None)
            curstr=mybam[mybam.iloc[:,0].str.contains('reads; of these:')].iloc[0,0]
            outstr=int(curstr.replace(" reads; of these:",""))
            bamsize[i]=outstr
        ## With bamfile size, normalize the read counts as a proportion
        repclass_norm=repclass_raw
        repname_norm=repname_raw
        repfamily_norm=repfamily_raw
        for i in sample_list:
            repclass_norm.loc[:,i]=repclass_raw[i]/bamsize[i]
            repname_norm.loc[:,i]=repname_raw[i]/bamsize[i]
            repfamily_norm.loc[:,i]=repfamily_raw[i]/bamsize[i]
        ## Save the normalized files
        repfamily_raw.to_csv("output/Repeat_expression_repfamily_raw.txt",sep="\t")
        repname_raw.to_csv("output/Repeat_expression_repname_raw.txt",sep="\t")
        repclass_raw.to_csv("output/Repeat_expression_repclass_raw.txt",sep="\t")
        repfamily_norm.to_csv("output/Repeat_expression_repfamily_norm.txt",sep="\t")
        repname_norm.to_csv("output/Repeat_expression_repname_norm.txt",sep="\t")
        repclass_norm.to_csv("output/Repeat_expression_repclass_norm.txt",sep="\t")


