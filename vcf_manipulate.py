import glob
import ntpath
import sys
import os
import subprocess

# Convert to annovar format
def convert_to_annovar(FilePath, patName, inFile, annovar):
    '''
    :param FilePath: Primary path of working directory
    :param patName: ID associated with a patient
    :param inFile: VCF file for neoantigens to be predicted from
    :param annovar: Dictionary housing annovar specific script locations. See README.md.
    :return: Name of the ANNOVAR ready vcf file for annovar_annnotation()
    '''
    print("INFO: Running convert2annovar.py on %s" % (inFile))

    # specify name for output file
    outDir = FilePath + "avready/"
    annovar_out_ready = outDir + patName + '.avinput'

    with open("logforannovarNeoPredPipe.txt", 'a') as logFile:
        cmd = ['perl', annovar['convert2annovar'], '-format', 'vcf4', inFile, '-outfile', annovar_out_ready, '-allsample',
             '-includeinfo', '-withfreq', '-comment']
        runconvert = subprocess.Popen(cmd, stdout=logFile, stderr=logFile)
        runconvert.wait()

    print('VCF Conversion Process complete %s'%(inFile))

    return(annovar_out_ready)

#Annotate with annovar
def annovar_annotation(FilePath, patName, inFile, annovar):
    '''
    :param FilePath: Primary path of working directory
    :param patName: ID associated with a patient
    :param inFile: ANNOVAR ready file.
    :param annovar: Dictionary housing annovar specific script locations. See README.md.
    :return: Name of ANNOVAR annotated exonic variant function file
    '''
    print("INFO: Running annotate_variation.pl on %s"%(inFile))

    # specify name for output file
    outDir = FilePath + "avannotated/"
    annovar_out_ready = outDir + patName + '.avannotated'

    with open("logforannovarNeoPredPipe.txt", 'a') as logFile:
        cmd = ['perl', annovar['annotatevariation'], '-out', annovar_out_ready, '-build', 'hg19', inFile, annovar['humandb'], '--comment']
        runannotate = subprocess.Popen(cmd, stdout=logFile, stderr=logFile)
        runannotate.wait()

    print('INFO: ANNOVAR annotation Process complete for %s' % (inFile))

    return(annovar_out_ready+".exonic_variant_function")

def get_coding_change(FilePath, patName, inFile, annovar):
    '''
    :param FilePath:
    :param patName:
    :param inFile:
    :param annovar:
    :return: Name of the fasta file with the coding change induced.
    '''
    print("INFO: Running coding_change.pl on %s" % (inFile))

    # specify name for output file
    outDir = FilePath + "fastaFiles/"
    coding_out_ready = outDir + patName + '.fasta'

    with open(coding_out_ready, 'a') as logfile:
        cmd = ['perl', annovar['coding_change'], inFile, annovar['gene_table'], annovar['gene_fasta'], '-includesnp', '-onlyAltering']
        runcoding = subprocess.Popen(cmd, stdout=logfile)
        runcoding.wait()

    print('INFO: Coding predictions complete for %s'%(inFile))
    return(coding_out_ready)

def predict_neoantigens(FilePath, patName, inFile, netMHCpan):
    '''
    :param FilePath:
    :param patName:
    :param inFile:
    :param netMHCpan:
    :return:
    '''
    pass