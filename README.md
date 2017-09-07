# ONT scripts
Scripts for population genetic analysis using Oxford Nanopore sequencing data.

## SNPs concordance analysis example
#### Calling SNPs via BCFtools
We first use SAMtools and BCFtools to call SNPs from Nanopore data. For Illumina reads, SNPs are called in the same way.

```
samtools mpileup -g -f Pvivax_Salvador_I.fasta sample1.ont.sorted.bam >sample1.ont.sorted.bcf
bcftools call -A -m sample1.ont.sorted.bcf >sample1.ont.sorted.Am.vcf
```

#### Quality filtering
SNPs are then filtered based on ```QUAL```. Here we only keep SNPs with ```QUAL``` greater than 20.

```
bcftools view -i "%QUAL>20" sample1.ont.sorted.Am.vcf -o sample1.ont.sorted.Am.q20.vcf
```

Indels are removed:

```
cat sample1.ont.sorted.Am.q20.vcf | grep -v "INDEL" >sample1.ont.sorted.Am.q20.snps.vcf
```

We use the ratio of ```DP4``` to ```DP``` as an additional filtering step. In this example, the ratio cutoff is 1. 

```
python filter_bcftools_dp4.py -i sample1.ont.sorted.Am.q20.snps.vcf -o sample1.ont.sorted.Am.q20.snps.dp4_1.vcf -r 1
```

#### SNP caller concordance
Now we can evaluate SNP concordance between Oxford Nanopore and Illumina runs. SNPs identified from Illumina reads are used as ground truth. 

```
python nanoporeSNPCaller_concordance.py -i sample1.iln.sorted.sorted.snps.vcf -n sample1.ont.sorted.Am.q20.snps.dp4_1.vcf -o sample1.sorted.Am.q20.snps.dp4_1.concor.txt -b -c bcftools
```

The output ```sample1.sorted.Am.q20.snps.dp4_1.concor.txt``` contains SNPs' position, read depth and binary classification information (TP/FP/TN/FN). Based upon these information, we are able to calculate precision and negative predictive value using the following command.

```
python nanopore_precision_npv.py -i sample1.sorted.Am.q20.snps.dp4_1.concor.txt -n sample1.npv.txt -p sample1.precision.txt -f sample1.freq.txt
```

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
We can then remove singletons. 
```
python filter_singleton_vcf.py -i merged_vcf.snps.miss40.vcf -o merged_vcf.snps.miss40.nosingletons.vcf
```

