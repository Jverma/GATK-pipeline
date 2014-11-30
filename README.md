GATK-pipeline
=============

A shell script which implements GATK pipeline for variant calling. Creates a <a href="https://www.broadinstitute.org/gatk/guide/article?id=4017"> GVCF</a> file for a sample BAM file. The GVCF files created by running the script on a cohort of samples can be jointly genotyped using <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php"> GenotypeGVCFs</a> 

<b>Download</b> GATK from https://www.broadinstitute.org/gatk/download

<b>Information</b> on GATK - https://www.broadinstitute.org/gatk/

<b>Dependencies</b>

* <a href="http://broadinstitute.github.io/picard/"> PicardTools</a>
* <a href="http://samtools.sourceforge.net/">SAMtools</a>

<b>Usage</b>

    sh GATK.sh
    
<b>Output</b>

    file_name_output.raw.snps.indels.g.vcf
    
    
GATK.pdf contains the explanation of the mathematics behind GATK-HaplotypeCaller.     
    





