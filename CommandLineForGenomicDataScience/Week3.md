# Answering Questions for Week 3  
## Loading programs  
module load samtools/1.2  
module load bcftools/1.2  
module load bowtie2/2.2.4  

### Questions  

#### Apply to questions 1 - 5: 
Generate a bowtie2 index of the wu_0_A genome using bowtie2-build, with the prefix 'wu_0'.
```
bowtie2-build wu_0.v7.fas wu_0
```

1.  How many sequences were in the genome?  
`grep ">" wu_0.v7.fas |wc -l`  
2. What was the name of the third sequence in the genome file? Give the name only, without the “>” sign.  
`grep ">" wu_0.v7.fas | head`  
3.  What was the name of the last sequence in the genome file? Give the name only, without the “>” sign.  
`grep ">" wu_0.v7.fas | tail`  
4.  How many index files did the operation create?  
`ls *.bt2 | wc -l`    
5. What is the 3-character extension for the index files created?

#### Apply to questions 6 - 14:   

Run bowtie2 to align the reads to the genome, under two scenarios: first, to report only full-length matches of the reads; and second, to allow partial (local) matches. 
All other parameters are as set by default. 

```
bowtie2 --local -x wu_0 -U wu_0_A_wgs.fastq -S Wu_0_A.local.sam
bowtie2 --end-to-end -x wu_0 -U wu_0_A_wgs.fastq -S Wu_0_A.E2E.sam
```
6.  How many reads were in the original fastq file?  
`echo $(cat *.fastq | wc -l)/4 | bc `

7.  How many matches (alignments) were reported for the original (full-match) setting? Exclude lines in the file containing unmapped reads.  
```
head bowtie_E2E.txt

147354 reads; of these:
  147354 (100.00%) were unpaired; of these:
    9635 (6.54%) aligned 0 times
    93780 (63.64%) aligned exactly 1 time
    43939 (29.82%) aligned >1 times
93.46% overall alignment rate

echo 93780+43939 | bc
```

8.  How many matches (alignments) were reported with the local-match setting? Exclude lines in the file containing unmapped reads.
```
head bowtie_local.txt

147354 reads; of these:
  147354 (100.00%) were unpaired; of these:
    6310 (4.28%) aligned 0 times
    84939 (57.64%) aligned exactly 1 time
    56105 (38.07%) aligned >1 times
    
echo 84939+56105 | bc
```

9.  How many reads were mapped in the scenario in Question 7?  
10.   How many reads were mapped in the scenario in Question 8?  

For the following set of questions (11 - 20), use the set of full-length alignments calculated under scenario 1 only. 
Convert this SAM file to BAM, then sort the resulting BAM file. 

```
samtools view -bS Wu_0_A.E2E.sam > Wu_0_A.E2E.bam
samtools sort Wu_0_A.E2E.bam -o Wu_0_A.sorted.E2E.bam
samtools index W_0_A.sorted.E2E.bam
```
11. How many reads had multiple matches in the scenario in Question 7? You can find this in the bowtie2 summary; note that by default bowtie2 only reports the best match for each read.  
12.  How many reads had multiple matches in the scenario in Question 8? Use the format above. You can find this in the bowtie2 summary; note that by default bowtie2 only reports the best match for each read.  
13.  How many alignments contained insertions and/or deletions, in the scenario in Question 7?  
`cat Wu_0_A.E2E.sam | cut -f 6 | egrep -c '[D]|[I]'`

14.  How many alignments contained insertions and/or deletions, in the scenario in Question 8?
`cat Wu_0_A.local.sam | cut -f 6 | egrep -c '[D]|[I]'`

#### Apply to questions 15 - 19:  

Compile candidate sites of variation using SAMtools mpileup for further evaluation with BCFtools. 
Provide the reference fasta genome and use the option “-uv” to generate the output in uncompressed VCF format for easy examination. 

```
samtools mpileup -f wu_0.v7.fas -uv Wu_0_A.sorted.E2E.bam > Wu_0_A.vcf
```
15.  How many entries were reported for Chr3?  
`cat *.vcf | grep -v "^#" | cut -f1 | grep -c "Chr3"`  
16.  How many entries have ‘A’ as the corresponding genome letter?  
`cat *.vcf | cut -f4 | grep -c "^A$"`  
17.  How many entries have exactly 20 supporting reads (read depth)?  
`cat *.vcf | cut -f8 | grep -c "DP=20"`  
18.  How many entries represent indels? 
`cat *.vcf | cut -f8 | grep "INDEL" | wc -l`  
19.  How many entries are reported for position 175672 on Chr1?  
`cat *.vcf | grep -v "^#" | grep "Chr1" | cut -f2  | grep -c "^175672$"` 


#### Apply to questions 20 - 24:   

Call variants using ‘BCFtools call’ with the multiallelic-caller model. For this, you will need to first re-run SAMtools mpileup with the BCF output format option (‘-g’)
to generate the set of candidate sites to be evaluated by BCFtools. 
In the output to BCFtools, select to show only the variant sites, in uncompressed VCF format for easy examination.

```
samtools mpileup -f wu_0.v7.fas -g Wu_0_A.sorted.E2E.bam > Wu_0_A.bcf
bcftools call -m -v -O v Wu_0_A.bcf > Wu_0_A.E2E.bcf.vcf
```
20.  How many variants are called on Chr3?  
`cat Wu_0_A.E2E.bcf.vcf | grep -v "^#" | cut -f1 | grep -c "Chr3"`  
21.  How many variants represent an A->T SNP? If useful, you can use ‘grep –P’ to allow tabular spaces in the search term.  
`cat Wu_0_A.E2E.bcf.vcf | grep -v "^#" | cut -f 4,5 | grep -P -c "^A\tT$"`  
22.  How many entries are indels?  
`cat Wu_0_A.E2E.bcf.vcf | grep -v "^#" | cut -f8 | grep "INDEL" | wc -l`  
23.  How many entries have precisely 20 supporting reads (read depth)?  
`cat Wu_0_A.E2E.bcf.vcf | grep -v "^#" | cut -f8 | grep "DP=20" | wc -l`  
24.  What type of variant (i.e., SNP or INDEL) is called at position 11937923 on Chr3?
`cat Wu_0_A.E2E.bcf.vcf | grep -v "^#" | grep -P "Chr3\t11937923"`  

