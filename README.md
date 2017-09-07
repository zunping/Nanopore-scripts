# ONT scripts
Scripts for population genetic analysis using Oxford Nanopore sequencing data.

## SNPs concordance analysis example
#### Calling SNPs via BCFtools
We first use SAMtools and BCFtools to call SNPs from Nanopore data. For Illumina reads, SNPs are called in the same way.

```samtools mpileup -g -f Pvivax_Salvador_I.fasta sample1.ont.sorted.bam >sample1.ont.sorted.bcf```<br />
```bcftools call -A -m sample1.ont.sorted.bcf >sample1.ont.sorted.Am.vcf```

#### Quality filtering
SNPs are then filtered based on ```QUAL``` and ```DP4``` values. Indels are removed.

```bcftools view -i "%QUAL>20" sample1.ont.sorted.Am.vcf -o sample1.ont.sorted.Am.q20.vcf```<br />
```cat sample1.ont.sorted.Am.q20.vcf | grep -v "INDEL" >sample1.ont.sorted.Am.q20.snps.vcf```
```python filter_bcftools_dp4.py -i sample1.ont.sorted.Am.q20.snps.vcf -o sample1.ont.sorted.Am.q20.snps.dp4_1.vcf```

#### SNP caller concordance
We now can evaluate SNP concordance between Oxford Nanopore and Illumina runs. SNPs identified from Illumina reads are used as ground truth. 

```python nanoporeSNPCaller_concordance.py -i sample1.iln.sorted.sorted.snps.vcf -n sample1.ont.sorted.Am.q20.snps.dp4_1.vcf -o sample1.sorted.Am.q20.snps.dp4_1.concor.txt -b -c bcftools```

## Variant calling analysis example
#### Merge multiple vcf files
First we combine several vcf files into one file by using GATK.

```java -jar GenomeAnalysisTK.jar -R 3D7_v25.fasta -T CombineVariants  --variant sample1.vcf --variant sample2.vcf --variant sample3.vcf --variant sample4.vcf --variant sample5.vcf -o merged_vcf.vcf -genotypeMergeOptions REQUIRE_UNIQUE ```

#### Remove INDELs and ReferenceInAll
Then we remove positions with `ReferenceInAll` or `INDEL`:
```cat merged_vcf.vcf | grep -Ev "ReferenceInAll|INDEL" >merged_vcf.snps.vcf```

#### Missingness filter
To apply the missingness filter, we run the following command:

```python filter_missingness_vcf.py -i merged_vcf.snps.vcf -o merged_vcf.snps.miss40.vcf -m 0.4```

This command will remove the SNPs with more than 0.4 missingness value and save the output to `merged_vcf.snps.miss40.vcf` . If there are five samples in total and three samples have `./.` for a certain SNP, then the missingness value is 0.6 and this SNP will be removed from the final output file. 

#### Filter singletons
We can then remove singletons. ```python filter_singleton_vcf.py -i merged_vcf.snps.miss40.vcf -o merged_vcf.snps.miss40.nosingletons.vcf```

