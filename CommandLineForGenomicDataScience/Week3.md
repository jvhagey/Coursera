# Answering Questions for Week 3
## Loading programs
module load samtools/1.2
module load bcftools/1.2
module load bowtie2/2.2.4

### Questions

Apply to questions 1 - 5: 
Generate a bowtie2 index of the wu_0_A genome using bowtie2-build, with the prefix 'wu_0'.
```
bowtie2-build wu_0.v7.fas wu_0
```

1.  
2.  
3.  
4.  
5.  

Apply to questions 6 - 14: 

Run bowtie2 to align the reads to the genome, under two scenarios: first, to report only full-length matches of the reads; and second, to allow partial (local) matches. 
All other parameters are as set by default. 

```
bowtie2 --local -x wu_0 -U wu_0_A_wgs.fastq -S Wu_0_A.local.sam
bowtie2 --end-to-end -x wu_0 -U wu_0_A_wgs.fastq -S Wu_0_A.E2E.sam
```

For the following set of questions (11 - 20), use the set of full-length alignments calculated under scenario 1 only. 
Convert this SAM file to BAM, then sort the resulting BAM file. 

```
samtools view -bS Wu_0_A.E2E.sam > Wu_0_A.E2E.bam
samtools sort Wu_0_A.E2E.bam -o Wu_0_A.sorted.E2E.bam
samtools index W_0_A.sorted.E2E.bam
```

Apply to questions 15 - 19: 

Compile candidate sites of variation using SAMtools mpileup for further evaluation with BCFtools. 
Provide the reference fasta genome and use the option “-uv” to generate the output in uncompressed VCF format for easy examination. 

```
samtools mpileup -uv Wu_0_A.sorted.E2E.bam > Wu_0_A.vcf
```

Apply to questions 20 - 24: 

Call variants using ‘BCFtools call’ with the multiallelic-caller model. For this, you will need to first re-run SAMtools mpileup with the BCF output format option (‘-g’)
to generate the set of candidate sites to be evaluated by BCFtools. 
In the output to BCFtools, select to show only the variant sites, in uncompressed VCF format for easy examination.

```
samtools mpileup -g Wu_0_A.sorted.E2E.bam > Wu_0_A.bcf
```
