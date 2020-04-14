
import sys
import re
import yaml

import os.path
import pandas as pd
# validate is only available from sbakemake 5.1
#from snakemake.utils import validate


def message(mes):
    sys.stderr.write("|--- " + mes + "\n")

def get_samples(sample_file):
    with open(sample_file) as fin:
        samples = [line.rstrip() for line in fin]
    return samples


def get_pe_chimeric_junction(wildcards):
    return os.path.join(config["staroutdir"], wildcards.sample, "pe", "Chimeric.out.junction")


def get_se_chimeric_junctions(wildcards):
    R1 = os.path.join(config["staroutdir"], wildcards.sample, "se", "R1", "Chimeric.out.junction")
    R2 = os.path.join(config["staroutdir"], wildcards.sample, "se", "R2", "Chimeric.out.junction")
    return {"R1": R1, "R2": R2}


def get_fastq(wildcards):
    fq1 = os.path.join(config['fastqdir'], "%s_1.fastq.gz" % wildcards.sample)
    fq2 = os.path.join(config['fastqdir'], "%s_2.fastq.gz" % wildcards.sample)
    return  {"fq1": fq1, "fq2": fq2}


configfile: "config.yaml"

samples = get_samples(config['samples'])

workdir: config['wdir']
message("The current working directory is " + config['wdir'])


rule all:
    input:
        expand("{sample}/circexplorer_circ.txt", sample=samples),
        expand("{sample}/circminer.circ_report", sample=samples),
        expand("{sample}/circexplorer2_circrna.txt", sample=samples),
        expand("{sample}/dcc/CircRNACount", sample=samples)


rule circexplorer:
    input:
        chimeric=get_pe_chimeric_junction,
        ref=config["ref"]["fasta"],
        genepred=config["ref"]["genepred"]
    output:
        circrna="{sample}/circexplorer_circ.txt",
        junctions=temp("{sample}_junction.txt")
    params:
        prefix="{sample}/circexplorer"
    log:
        parse_stderr="logs/{sample}.circexplorer_parse.e",
        parse_stdout="logs/{sample}.circexplorer_parse.o",
        explore_stderr="logs/{sample}.circexplorer_explore.e",
        explore_stdout="logs/{sample}.circexplorer_explore.o",
    shell:
        """
        set +u; \
        source /home/faraut/dynawork/CircRNA/softwares/CIRCexplorer/circenv/bin/activate ;
        set -u
        circeplorerpath=/home/faraut/dynawork/CircRNA/softwares/CIRCexplorer/circ
        python $circeplorerpath/star_parse.py {input.chimeric} \
               {output.junctions} 1>{log.parse_stdout} &>{log.parse_stderr}
        python $circeplorerpath/CIRCexplorer.py -j {output.junctions}  \
               -g {input.ref} -r {input.genepred} -o {params.prefix} \
                1>{log.explore_stdout} &>{log.explore_stderr}
        """


rule circexplorer2:
    input:
        unpack(get_fastq),
        ref=config["ref"]["fasta"],
        gtf=config["ref"]["gtf"],
        genepred=config["ref"]["genepred"]
    output:
        circrna="{sample}/circexplorer2_circrna.txt",
        junctions="{sample}/circexplorer2_BSJ.bed"
    log:
        align="logs/{sample}.circexplorer2_align.log",
        parse="logs/{sample}.circexplorer2_parse.log",
        annotate="logs/{sample}.circexplorer2_annotate.log"
    params:
        bowtie1_index = config['bowtie1_index'],
        bowtie2_index = config['bowtie2_index'],
    shell:
        """
        set +u; \
        source /home/faraut/dynawork/CircRNA/softwares/CIRCexplorer2/circ2env/bin/activate \
        set -u
        CIRCexplorer2 align  -G {input.gtf} -i {params.bowtie1_index} \
                      -j {params.bowtie2_index} -f {input.fq1},{input.fq2} \
                      -b {output.junctions} \
                      > {log.align}
        CIRCexplorer2 annotate -r {input.genepred} -g {input.ref}  \
           -b {output.junctions} -o {output.circrna} > {log.annotate}
        """
#CIRCexplorer2 parse -t STAR {input.chimeric} -b {output.junctions}> {log.parse} > {log.parse}

rule circminer:
    input:
        unpack(get_fastq),
        ref=config["ref"]["fasta"],
        gtf=config["ref"]["gtf"]
    output:
        "{sample}/circminer.circ_report",
        "{sample}/circminer.candidates.pam"
    params:
        prefix="{sample}/circminer"
    log:
        stderr="logs/{sample}.circminer.e",
        stdout="logs/{sample}.circminer.o"
    shadow:
        "shallow"
    shell:
        """
        module load compiler/gcc-5.3.0
        circminer=/work2/genphyse/dynagen/tfaraut/CircRNA/softwares/circminer/circminer
        $circminer -r {input.ref} -g {input.gtf} -s {input.fq1} -2 {input.fq2} \
        -o {params.prefix} -k 20 -a 0 -S 500 1>{log.stdout} 2>{log.stderr}
        """

rule dcc:
    input:
        unpack(get_se_chimeric_junctions),
        pechim=get_pe_chimeric_junction,
        ref=config["ref"]["fasta"],
        gtf=config["ref"]["gtf"]
    output:
        "{sample}/dcc/CircRNACount",
        "{sample}/dcc/CircCoordinates",
        "{sample}/dcc/CircSkipJunctions"
    log:
        stderr="logs/{sample}.dcc.e",
        stdout="logs/{sample}.dcc.o"
    threads:
        8
    shadow: "shallow"
    params:
        temp="{sample}/dcc/DCC_TMP",
        outputdir="{sample}/dcc"
    shell:
        """
        set +u; \
        source /work2/genphyse/dynagen/tfaraut/CircRNA/softwares/DCC/DCC/dccenv/bin/activate ; \
        set -u
        DCC -t {params.temp} -D -mt1 {input.R1} -mt2 {input.R2} -T {threads} \
            -an {input.gtf} -Pi -F -N -Nr 1 1 -fg -A {input.ref} {input.pechim} \
            -O {params.outputdir} 1>{log.stdout} 2>{log.stderr}
        """
