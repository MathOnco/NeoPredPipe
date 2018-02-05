#!/usr/bin/env python2.7

#runs the proper epitope predictions for topiary for MHCI molecules

import glob
import ntpath
import sys
import os
import subprocess

#input folder containing *.vcf files
arg1 = sys.argv[1]
#output folder new_head
arg2 = sys.argv[2]
#alleles folder
#arg3 = sys.argv[3]

print('Begin...')

#get list of vcf files and complete directories
def get_input_vcf(arg1):
	wild = '*.txt'
	vcf = arg1 + wild

	#generate list of files for inputs
	vcf_in = glob.glob(vcf)

	vcf_out = dict()
	for item in vcf_in:
		item2 = item
		vcf_out.update({item:item2})
	print('Input files found.')
	return vcf_out

#requires input, output, and the file list generated from get_input_vcf
#adds appropriate headers to files
def update_headers(arg1,arg2, inputs):
	#make sure all files have proper vcf headers
	for item in inputs:

		ind = item.rfind('/')
		outname1 = item[ind+1:len(item)].split('.vcf')
		outname =  os.getcwd() + '/' + arg2 + outname1[0] + '_new_head.vcf'

		with open(inputs[item], 'r') as vcfin:
			with open(outname, 'w') as vcfout:
				vcfout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + str(outname1[0]) + '\n')

				for line in vcfin:
					if line.startswith('##INFO='):
						vcfout.write(line)
						vcfout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + str(outname1[0]) + '\n')
					else:
						vcfout.write(line)
	print('Headers have been updated.')

#Convert the new header files to annovar format
def convert_to_annovar(arg2):
	print 'Converting to ANNOVAR format'
	with open('logfile.txt', 'a') as logfile:
		#Gets absolute path to the new header files and annovar
		arg2_full = os.path.abspath(arg2) + '/'

		#gets all the files into a list that need converting to annovar format
		for_convert = glob.glob(arg2_full + '*.vcf')

		#get the basename for the annovar converted files
		for_convert_dict = dict()
		for item in for_convert:
			file_name = item[item.rfind('/')+1:item.rfind('.')]
			for_convert_dict.update({file_name:item})

		annovar_full = '/Users/Ryan/Desktop/Bioinformatics_Tools/annovar/convert2annovar.pl'
		annovar_out = '/Users/Ryan/Desktop/TEST/av_ready/'

		for item in for_convert_dict:
			print 'Running convert2annover.pl: ' + item + '\n'

			#specify name for output file
			annovar_out_ready = annovar_out + item + '.avinput'

			#runconvert = subprocess.Popen(['perl', annovar_full, '-format', 'vcf4', for_convert_dict[item], '-outfile', annovar_out_ready,
			#	'-allsample', '-includeinfo', '-withfreq','-comment'], shell=True)
			#with open('logfile.txt', 'a') as logfile:
			runconvert = subprocess.Popen(['perl', annovar_full, '-format', 'vcf4', for_convert_dict[item], '-outfile', annovar_out_ready, '-allsample', '-includeinfo', '-withfreq','-comment'], stdout=logfile)
			runconvert.wait()
			print '\nVCF Conversion Process complete\n'
	return annovar_out


#Annotate with annovar
def annovar_annotation(av_ready):
	print 'Annotating with annovar'
	with open('logfile.txt', 'a') as logfile:
		#output directory for final annovar annotated files
		annotated_out = '/Users/Ryan/Desktop/TEST/av_annotated/'
		annovar_annotate_script = '/Users/Ryan/Desktop/Bioinformatics_Tools/annovar/annotate_variation.pl'

		#get filenames of the new_head.vcf files in av_ready folder
		for_annotation = glob.glob(av_ready + '*.avinput')

		#get the basename for the annovar converted files
		#key is stripped name
		#value is the file with file path
		for_annotation_dict = dict()
		for item in for_annotation:
			file_name = item[item.rfind('/')+1:item.rfind('.')]
			for_annotation_dict.update({file_name:item})

		#perform perl annovar annotation script
		for item in for_annotation_dict:
			print 'Running annotate_variation.pl:' + item +'\n'

			#prepare output
			annovar_out_ready = annotated_out + item

			runannotate = subprocess.Popen(['perl',annovar_annotate_script, '-out', annovar_out_ready, '-build', 'hg19', for_annotation_dict[item],'/Users/Ryan/Desktop/Bioinformatics_Tools/annovar/humandb/', '--comment'], stdout = logfile)
			runannotate.wait()

			print '\nANNOVAR annotation Process complete\n'



#files to predict epitopes from
inputs = get_input_vcf(arg1)

#update headers to proper vcf headers
update_headers(arg1,arg2,inputs)

#converts the new header files to annovar input format files, yields a path to the input for annovar_annotation
annovar_out = convert_to_annovar(arg2)

#Uncomment if the process above has been completed
#annovar_out = '/Users/Ryan/Desktop/TEST/av_ready/'

#annotate with annovar (annovar_out is from conversion)
annovar_annotation(annovar_out)


'''
#USE THIS FOR REFSEQ ANNOTATION
perl annotate_variation.pl -out test.txt -build hg19 example/multiSNV_4_filtered_annotated_new_head.avinput humandb/ --comment
'''
'''
#USE THIS FOR ENSEMBL TRANSCRIPTS
perl annotate_variation.pl -out test.txt -build hg19 example/multiSNV_4_filtered_annotated_new_head.avinput humandb/ --comment -dbtype ensGene
'''

