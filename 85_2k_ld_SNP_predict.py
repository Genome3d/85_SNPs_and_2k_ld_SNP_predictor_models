#!/usr/bin/env python Daniel Ho 08/05/2022

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import pickle


#Input file is assumed to be plink raw format genereated by "plink --recode A" with plink format files (bed, bim and fam).
#Please refer to https://www.cog-genomics.org/plink/1.9/formats#raw for the raw format information.
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("file", help="enter the plink raw genotype file name")
args = parser.parse_args()

fullfilename = args.file
pos = fullfilename.find('.raw')
filename = fullfilename[:pos]


# load 85QTLs_2k_ld.snps
def load_table(datafile):

	SNP_alleles = []

	inputfile = open(datafile, 'r')
	while True:
		input_line = inputfile.readline()
		input_line = input_line[:-1]
		if not input_line:
			break
		else:
			SNP_alleles.append(input_line)
	inputfile.close()
	return SNP_alleles



#check the SNP-allele specification fo the input plink raw file with 85QTLs_2k_ld.snps.
def Correct_SNP_alleles(datafile,SNP_alleles):

	Correct = False
	SNP_table = pd.read_csv(datafile,sep=" ")
	x_features = SNP_table[SNP_table.columns[6:]]
	if len(x_features.columns) == len(SNP_alleles):
		for x in x_features.columns:
			if x not in SNP_alleles:
				print("Error: incorrect SNP allele " + x)
				break
		Correct = True
	else:
		print ("Error: incorrect SNP number")
	return Correct




#create two tissue specific eQTL tables from the input raw files
def Create_eQTL_matrix(datafile):
	eQTL_table = pd.read_csv(datafile,sep=" ")
	for x in [1,2]:
		if x == 1:
			TSG_list = pd.read_csv('data/Fetal_Tissue_SNP_Gene_sig_mapping_sorted_beta.txt',sep='\t')
		else:
			TSG_list = pd.read_csv('data/Colon_Tissue_SNP_Gene_sig_mapping_sorted_beta.txt',sep='\t')
		tmp_eQTL_table = eQTL_table.copy()
		samples_header = list (eQTL_table.columns[0:6])
		SNP_w_list = list(eQTL_table.columns[6:])
		tmp_SNP_w_list = SNP_w_list.copy()
		TSG_cols = []
		drop_SNPs = []
		for rsSNP_w in tmp_SNP_w_list :
			rsSNP = rsSNP_w.split('_')[0]
			tmp_TSGs = list ( TSG_list.Tissue[TSG_list.SNP == rsSNP] + "--" + rsSNP_w + "--" + TSG_list.Gene_Name[TSG_list.SNP == rsSNP])
			tmp_Effects = list ( TSG_list.Effect_Size[TSG_list.SNP == rsSNP])
			if len(tmp_Effects) != 0 :
				count = 0
				for TSG in tmp_TSGs:
					tmp_eQTL_table[TSG] = tmp_eQTL_table[rsSNP_w] * tmp_Effects[count]
					count = count + 1
					TSG_cols.append(TSG)
				drop_SNPs.append(rsSNP_w)
		for rsSNP_w in drop_SNPs:
			SNP_w_list.remove(rsSNP_w)
		TSG_col_sorted = sorted(TSG_cols)
		new_eQTL_table = tmp_eQTL_table[samples_header + TSG_col_sorted]
		if x == 1:
			new_eQTL_table.to_csv('data/Feta_Tissue_SNP_Gene_eQTL_table_beta.txt', sep='\t', index=False)
		else:
			new_eQTL_table.to_csv('data/Colon_Tissue_SNP_Gene_eQTL_table_beta.txt', sep='\t', index=False)
			
	





#load three predictor models from the data directory to generate the AUC prediction results
def Prediction_results(datafile):
	filename = 'data/85_2k_LD_SNP_predictor2.sav'
	eQTL_table = pd.read_csv(datafile,sep=" ")
	X = eQTL_table[eQTL_table.columns[6:]]
	Y = eQTL_table['PHENOTYPE'] - 1
	loaded_model = pickle.load(open(filename, 'rb'))
	result = roc_auc_score(Y, loaded_model.predict_proba(X)[:,1])
	print("AUC prediction result for SNP only model : " + "{:.2f}".format(result))
    
	filename = 'data/85_2k_LD_SNP_adult_fetal_brain_predictor2.sav'
	eQTL_table = pd.read_csv("data/Feta_Tissue_SNP_Gene_eQTL_table_beta.txt",sep="\t")
	X = eQTL_table[eQTL_table.columns[6:]]
	Y = eQTL_table['PHENOTYPE'] - 1
	loaded_model = pickle.load(open(filename, 'rb'))
	result = roc_auc_score(Y, loaded_model.predict_proba(X)[:,1])
	print("AUC prediction result adult and fetal brain tissue model : " + "{:.2f}".format(result))

	filename = 'data/85_2k_LD_SNP_colon_predictor2.sav'
	eQTL_table = pd.read_csv("data/Colon_Tissue_SNP_Gene_eQTL_table_beta.txt",sep="\t")
	X = eQTL_table[eQTL_table.columns[6:]]
	Y = eQTL_table['PHENOTYPE'] - 1
	loaded_model = pickle.load(open(filename, 'rb'))
	result = roc_auc_score(Y, loaded_model.predict_proba(X)[:,1])
	print("AUC prediction result colon tissue model : " + "{:.2f}".format(result))



interlines = []

SNP_alleles = load_table("data/85QTLs_2k_ld.snps")

if Correct_SNP_alleles(fullfilename,SNP_alleles):
	Create_eQTL_matrix(fullfilename)
	Prediction_results(fullfilename)

