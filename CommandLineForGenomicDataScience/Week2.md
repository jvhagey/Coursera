# Answering question for week 2

## Loading programs
module load samtools/1.2  
module load BEDTools/2.23.0  

## Running 
samtools flagstat athal_wu_0_A.sorted.bam  #Q1  

samtools sort athal_wu_0_A.bam athal_wu_0_A.sorted  
samtools index athal_wu_0_A.sorted.bam  

1.  How many alignments does the set contain?  
`samtools flagstat athal_wu_0_A.bam`

2.  How many alignments show the read’s mate unmapped?  
`samtools view athal_wu_0_A.bam | cut -f 7 | grep -c '*'`

3.  How many alignments contain a deletion (D)?  
`samtools view athal_wu_0_A.bam | cut -f 6 | grep -c "D"`

4.  How many alignments show the read’s mate mapped to the same chromosome?  
`samtools view athal_wu_0_A.bam | cut -f 7 | grep -c '='`

5.  How many alignments are spliced?  
`samtools view athal_wu_0_A.bam | cut -f 6 | grep -c "N"`

6.  How many alignments does the set contain?  
`samtools view -H athal_wu_0_A.bam | grep "SN"`


7. How many alignments does the set contain?  
```
samtools flagstat extracted.sam 
samtools index athal_wu_0_A.bam
samtools view -h athal_wu_0_A.bam Chr3:11777000-11794000 > extracted.sam
```

8. How many alignments show the read’s mate unmapped?  
`samtools view extracted.bam | cut -f 7 | grep -c '*'`

9.  How many alignments contain a deletion (D)?  
`samtools view extracted.bam | cut -f 6 | grep -c 'D'`

How many alignments show the read’s mate mapped to the same chromosome?  
`samtools view athal_wu_0_A.bam | cut -f 7 | grep -c '='`

10.  How many alignments are spliced?  
`samtools view athal_wu_0_A.bam | cut -f 7 | grep -c 'N'`

11.  How many sequences are in the genome file?  
`samtools view -H athal_wu_0_A.bam | grep -c "SN:"`

12.  What is the length of the first sequence in the genome file?  
`samtools view -H athal_wu_0_A.bam | grep "SN:"`

13.  What alignment tool was used?  
`samtools view -H athal_wu_0_A.bam | grep "^@PG"`

14.  What is the read identifier (name) for the first alignment?    
`samtools view athal_wu_0_A.bam | cut -f 1 | head`

15.  What is the start position of this read’s mate on the genome?  
```
#Give this as ‘chrom:pos’ if the read was mapped, or ‘*” if unmapped.
samtools view athal_wu_0_A.bam | cut -f 3,8 | head
```

16.  How many overlaps (each overlap is reported on one line) are reported?   
```
bedtools bamtobed -i athal_wu_0_A.bam > data.bed
bedtools intersect -wo -a athal_wu_0_A_annot.gtf -b data.bed | wc -l 
OR
bedtools intersect -wo -a athal_wu_0_A_annot.gtf -b athal_wu_0_A.bam | wc -l
```

17.  How many of these are 10 bases or longer?  
```
bedtools intersect -wo -a athal_wu_0_A_annot.gtf -b athal_wu_0_A.bed | cut -f16 > count.txt
awk '$1>9{c++} END{print c+0}' count.txt
```


18.  How many alignments overlap the annotations?  
`bedtools intersect -wo -a athal_wu_0_A_annot.gtf -b data.bed | wc -l`

19.  Conversely, how many exons have reads mapped to them?  
```
bedtools intersect -wo -a athal_wu_0_A_annot.gtf -b data.bed | cut -f4 | sort -u | wc -l
OR
bedtools intersect -wo -a athal_wu_0_A_annot.gtf -b data.bed | cut -f5 | sort -u | wc -l
```

20.  If you were to convert the transcript annotations in the file “athal_wu_0_A_annot.gtf” into BED format, how many BED records would be generated?  
`bedtools intersect -wo -a athal_wu_0_A_annot.gtf -b data.bed | cut -f9 | cut -d " " -f4 | sort -u | wc -l`
