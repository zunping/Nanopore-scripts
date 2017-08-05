# ONT scripts
Scripts for population genetic analysis using Oxford Nanopore sequencing data.
## Missingness filter
First we combined several vcf files into one file by using GATK.

```java -jar GenomeAnalysisTK.jar -R 3D7_v25.fasta -T CombineVariants  --variant sample1.vcf --variant sample2.vcf --variant sample3.vcf --variant sample4.vcf -o merged_vcf.vcf -genotypeMergeOptions REQUIRE_UNIQUE ```

To apply the missingness filter

```python filter_missingness_vcf.py -i merged_vcf.snps.vcf -o merged_vcf.snps.miss40.vcf -m 0.4```
