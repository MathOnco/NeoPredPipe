import glob
import ntpath
import sys
import os
import subprocess
import shutil
from Bio import SeqIO

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
    :param FilePath: Primary path of working directory
    :param patName: ID associated with a patient
    :param inFile: ANNOVAR exonic variant function file.
    :param annovar: Dictionary housing annovar specific script locations. See README.md.
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

def ReformatFasta(inFile):
    '''
    Reformats header for netMHCpan of all fasta files

    :param inFile: Fasta file created by coding_change.pl
    :return: Name of new fasta file.
    '''
    newFasta = inFile.rstrip(".fasta") + ".reformat.fasta"
    with open(inFile, 'r') as fasta:
        with open(newFasta, 'w') as outFasta:
            for line in fasta.readlines():
                if '>' in line:
                    outFasta.write(line.replace(" ",";"))
                else:
                    outFasta.write(line)

    return(newFasta)

def ExtractSeq(seq_record, pos, n):
    '''
    Extracts the proper range of amino acids from the sequence and the epitope length

    :param seq: Sequence of mutated amino acid sequence.
    :param pos: Location of altered amino acid.
    :param n: Length of epitope for predictions.
    :return: Returns an amino acid sequence of appropriate lengths.
    '''
    seq = str(seq_record.seq)
    if pos<n: # Cases where start is less than n
        miniseq = seq[0:pos+(n)]
    elif len(seq)-pos < n: # Cases where end is less than n
        miniseq = seq[pos - (n - 1):len(seq)]
    else:
        miniseq = seq[pos-(n-1):pos+(n)]
    return(miniseq)

def MakeTempFastas(inFile, epitopeLens):
    '''
    Creates a peptide sequence fasta for netMHCpan containing only the region of interest. No WILDTYPE, immediate-stopgain, or stop-loss.

    :param inFile: Fasta file with reformatted headers
    :param epitopeLens: User option of epitope lengths for predictions
    :return: Not sure yet.
    '''
    eps = {n: 0 for n in epitopeLens}
    for n in epitopeLens:
        mySeqs = []
        for seq_record in SeqIO.parse(inFile, 'fasta'):
            if 'wildtype' not in seq_record.id.lower() and 'immediate-stopgain' not in seq_record.id.lower() and 'from;*;to;' not in seq_record.id.lower():
                pos = int(seq_record.id.replace(";;",";").split(";")[5])-1

                miniseq = ExtractSeq(seq_record, pos, n)
                mySeqs.append(">"+seq_record.id+"\n"+miniseq+"\n")
        eps[n] = mySeqs

    tmpFiles = {}
    for n in epitopeLens:
        tmpFasta = inFile.rstrip('.reformat.fasta') + ".tmp.%s.fasta"%(n)
        tmpFiles.update({n:tmpFasta})
        with open(tmpFasta, 'w') as outFile:
            for line in eps[n]:
                outFile.write(line)

    return(tmpFiles)

def ConstructAlleleHelper(s):
    return(s[:4].lower() + s[4:].capitalize())

def ConstructAlleles(hlas):
    '''
    Constructs the proper HLA input from HLA calls.

    :param hlas: list of HLA types for the Patient
    :return: list of normalized HLA identifiers for netMHCpan
    '''
    with open("netMHCpanAlleles.txt",'r') as alleles:
        allAlleles = [i.rstrip('\n').lower() for i in alleles.readlines()]

    hlas = [i.replace("hla_","hla-") for i in hlas]
    hlas = [hla.replace("_","",1) for hla in hlas if 'NA' not in hla]
    hlas = [hla.replace("_",":",1) for hla in hlas if 'NA' not in hla]

    netMHCpanHLAS = []
    for hla in hlas:
        if hla.replace("_","")[0:-2] in allAlleles:
            netMHCpanHLAS.append(hla.replace("_","")[0:-2].upper())
        elif hla.replace("_","")[0:-4] in allAlleles:
            netMHCpanHLAS.append(hla.replace("_", "")[0:-4].upper())
        else:
            sys.exit("HLA type not found.")

    return(list(set(netMHCpanHLAS)))

def predict_neoantigens(FilePath, patName, inFile, hlas, epitopeLens, netMHCpan):
    '''
    Strips out all WILDTYPE and IMMEDIATE-STOPGAIN from fasta file.

    :param FilePath: Primary path of working directory
    :param patName: ID associated with a patient
    :param inFile: Fasta file with reformatted coding changes.
    :param hlas: HLA types for the patient.
    :param epitopeLens: List of epitope lengths to predict
    :param netMHCpan: Dictionary housing netMHCpan specific script locations and data. See README.md.
    :return: netMHCpan
    '''

    hlasnormed = ConstructAlleles(hlas)

    print("INFO: Predicting neoantigens for %s" % (patName))

    for n in epitopeLens:
        output_file = FilePath +'tmp/%s.epitopes.%s.txt' % (patName, n)
        with open(output_file, 'a') as epitope_pred:
            print("INFO: Running Epitope Predictions for %s on epitopes of length %s"%(patName,n))
            cmd = ['netMHCpan', '-l', str(n), '-a', ','.join(hlasnormed), '-f', inFile[n]]
            netMHC_run = subprocess.Popen(cmd, stdout=epitope_pred, stderr=epitope_pred)
            netMHC_run.wait()

    print("INFO: Predictions complete for %s on epitopes of length %s" % (patName, n))

