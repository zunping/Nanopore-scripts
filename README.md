# ONT scripts
Scripts for population genetic analysis using Oxford Nanopore sequencing data.

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

