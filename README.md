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
module load bioinfo/snakemake-4.8.0
```

```bash
snakemake --jobs 30 --cluster-config cluster.yaml --drmaa " --mem-per-cpu={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" -p -n
```
