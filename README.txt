# python environment:
# 85_2k_ld_SNP_predict.py was developed with python3.7.3 with modules: numpy, pandas, sklearn.metrics, pickle, and argparse (required to be installed in the
# python environment for 85_2k_ld_SNP_predict.py to work properly. 
# Example: python 85_2k_ld_SNP_predict.py GenotypeData.raw
# AUC prediction result for SNP only model : 0.90
# AUC prediction result adult and fetal brain tissue model : 0.45
# AUC prediction result colon tissue model : 0.70

# The PHENOTYPE is defined using plink format that "1" is control and "2" is case. 

# Input file of 85_2k_ld_SNP_predict.py is assumed to be plink raw format genereated by "plink --recode A"
# with plink format files (bed, bim and fam).
# Please refer to https://www.cog-genomics.org/plink/1.9/formats#raw for the raw format information.
# If the order of a1 and a2 of a SNP allele is incorrect, you could try to fix the order problem by "plink --a1-allele"

# The input raw format file is required to have 232 data columns of the SNPs with their alleles specified
# in 85QTLs_2k_ld.snps file. Please check 85QTLs_2k_ld.snps file for the SNP-allele specification.

# 1, 85_2k_ld_SNP_predict.py check the SNP-allele specification fo the input plink raw file with 85QTLs_2k_ld.snps.
# 2, 85_2k_ld_SNP_predict.py create two tissue specific eQTL tables from the input raw files 
#    using two tissue specific mapping files in the data directory.
# 3, 85_2k_ld_SNP_predict.py load three predictor models from the data directory to generate the AUC prediction results
#    for using SNPs only, adult and fetal eQTLs and Colon eQTLs.

# The original 10 fold AUC results as the following:
# for SNPs only: AUC 0.79
# for adult and fetal brain: AUC 0.69 
# for colon tissue: AUC 0.68

