# circbenchmark

```bash
module load bioinfo/samtools-1.9
module load bioinfo/bcftools-1.9
module load bioinfo/seqtk-1.2
module load bioinfo/bwa-0.7.17
module load bioinfo/STAR-2.6.0c
module load bioinfo/RSEM-1.3.3
module load bioinfo/bedops-v2.4.35
module load bioinfo/bowtie-1.2.1.1
module load bioinfo/bowtie2-2.3.5.1
module load bioinfo/tophat-2.1.1
module load bioinfo/snakemake-4.8.0
```

```bash
snakemake --jobs 30 --cluster-config cluster.yaml --drmaa " --mem-per-cpu={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" -p -n
```

#### Merging HeLa RNase R- and RNase R+ fastq files for circminer

RNase R-
```bash
FASTQDIR=/work2/genphyse/dynagen/tfaraut/CircRNA/benchmarking/CircMinerPaper/data/fastq

for strand in 1 2;
do
   for replicate in  SRR1637089 SRR1637090;
   do
     cat ${FASTQDIR}/${replicate}_${strand}.fastq.gz
   done > ${FASTQDIR}/HeLaRminus_${strand}.fastq.gz
   for replicate in SRR1636985 SRR1636986;
   do
     cat ${FASTQDIR}/${replicate}_${strand}.fastq.gz
   done > ${FASTQDIR}/HeLaRplus_${strand}.fastq.gz
done


```
