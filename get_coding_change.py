#!/usr/bin/env python2.7

#Get fasta files for ANNOVAR annotated scripts.

import sys
import subprocess
import glob
import os

#input files (ANNOVAR exonic_variant_function) directory
arg1 = sys.argv[1]
#output directory
arg2 = sys.argv[2]

#obtain the protein sequence for the snps
def get_coding_change(arg1, arg2):
	arg1_full = os.path.abspath(arg1) + '/'
	arg2_full = os.path.abspath(arg2) + '/'

	coding_full = '/Users/Ryan/Desktop/Bioinformatics_Tools/annovar/coding_change.pl'
	gene_table = '/Users/Ryan/Desktop/Bioinformatics_Tools/annovar/humandb/hg19_refGene.txt'
	gene_fasta = '/Users/Ryan/Desktop/Bioinformatics_Tools/annovar/humandb/hg19_refGeneMrna.fa'
	for_coding = glob.glob(arg1_full + '*.exonic_variant_function')


	#get the basename for the annovar converted files
	for_coding_dict = dict()
	for item in for_coding:
		file_name = item[item.rfind('/')+1:item.rfind('.')]
		for_coding_dict.update({file_name:item})

	for item in for_coding_dict:
		print 'Running coding_change.pl:' + item + '\n'

		output_name = arg2_full + item + '.fasta'
		with open(output_name, 'a') as logfile:
			runcoding = subprocess.Popen(['perl',coding_full, for_coding_dict[item], gene_table, gene_fasta,'-includesnp','-onlyAltering'], stdout=logfile)
			runcoding.wait()

		print 'Coding predictions complete'

get_coding_change(arg1,arg2)