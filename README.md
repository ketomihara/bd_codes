# bd_codes
Scripts for our paper (Tomihara & Kiuchi, 2023, bioRxiv. doi: 10.1101/2023.04.01.535244) are provided here.

**Mapping_VariantCalling.sh**  
Dependencies: STAR 2.7.20, GATK 4.3.0, samtools 1.6, sambamba 1.0.0. Other versions may be also fine.  
Genome assembly and gene annotation files of Bombyx mori can be obtained from SilkBase (https://silkbase.ab.a.u-tokyo.ac.jp/) or KAIKObase (https://kaikobase.dna.affrc.go.jp/).  
RNA-seq reads were provided from Wu et al. 2016, Sci. Rep. doi: 10.1038/srep26114.
This codes requires RNA_processing_list.txt, which specifies the location of the data files.  

**VCF_SNPs.R**  
Dependencies: R 3.6.3, dplyr, ggplot2, stringr, reshape2. Other versions may be also fine.  
This codes counts heterozygous SNPs from VCF file obtained from Mapping_VariantCalling.sh.
