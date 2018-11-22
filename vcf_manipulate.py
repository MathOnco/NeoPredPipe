#!/usr/bin/env python

'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
Contributions from: Eszter Lakatos
'''

import sys
import os
import subprocess
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
    outDir = "avready/"
    annovar_out_ready = outDir + patName + '.avinput'

    with open("logforannovarNeoPredPipe.txt", 'a') as logFile:
        cmd = ['perl', annovar['convert2annovar'], '-format', 'vcf4', inFile, '-outfile', annovar_out_ready, '-allsample',
             '-includeinfo', '-withfreq', '-comment']
        runconvert = subprocess.Popen(cmd, stdout=logFile, stderr=logFile)
        runconvert.wait()

    print('INFO: VCF Conversion Process complete %s'%(inFile))

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
    outDir = "avannotated/"
    annovar_out_ready = outDir + patName + '.avannotated'

    with open("logforannovarNeoPredPipe.txt", 'a') as logFile:
        cmd = ['perl', annovar['annotatevariation'], '-out', annovar_out_ready, '-build', annovar['build'], inFile, annovar['humandb'], '--comment']
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
    outDir = "fastaFiles/"
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
    newFasta = inFile.replace(".fasta",".reformat.fasta")
    with open(inFile, 'r') as fasta:
        with open(newFasta, 'w') as outFasta:
            for line in fasta.readlines():
                if '>' in line:
                    outFasta.write(line.replace(" ",";"))
                else:
                    outFasta.write(line)

    return(newFasta)

def ExtractSeq(seq_record, pos, n, frameshift=False):
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

    if frameshift:
        if pos<n:
            miniseq = seq[0:len(seq)] # When start is not n away from pos, we go until the end of the sequence
        else:
            miniseq = seq[pos-(n-1):len(seq)] # All other cases, we still go until the end of the sequence

    return(miniseq)

def MakeTempFastas(inFile, epitopeLens):
    '''
    Creates a peptide sequence fasta for netMHCpan containing only the region of interest. No WILDTYPE, immediate-stopgain, or stop-loss.

    :param inFile: Fasta file with reformatted headers
    :param epitopeLens: User option of epitope lengths for predictions
    :return: Not sure yet.
    '''
    eps = {n: 0 for n in epitopeLens}
    epsIndels = {n: 0 for n in epitopeLens}
    for n in epitopeLens:
        mySeqs = []
        mySeqsIndels = []
        for seq_record in SeqIO.parse(inFile, 'fasta'):
            if 'wildtype' not in seq_record.id.lower() and 'immediate-stopgain' not in seq_record.id.lower() and 'from;*;to;' not in seq_record.id.lower() and 'silent' not in seq_record.id.lower() and 'fs*' not in seq_record.id.lower() and 'delins' not in seq_record.id.lower():
                # TODO: Add a regex expression to extract the position since it's variable with versions of ANNOVAR
                # TODO: Regex code for this might be r"\w*((?i)position;\d+;(?-i))\W*"
                try:
                    pos = int(seq_record.id.replace(";;",";").split(";")[5])-1
                except ValueError:
                    pos = int(seq_record.id.replace(";;",";").split(";")[6])-1

                if 'dup' in seq_record.id.lower() or 'del' in seq_record.id.lower() or 'ins' in seq_record.id.lower(): #if mutation is potentially frameshift
                    miniseq = ExtractSeq(seq_record, pos, n, True) # flag for frameshift
                    mySeqsIndels.append(">"+seq_record.id+"\n"+miniseq+"\n")
                else:
                    miniseq = ExtractSeq(seq_record, pos, n)
                    mySeqs.append(">"+seq_record.id+"\n"+miniseq+"\n")

        eps[n] = mySeqs
        epsIndels[n] = mySeqsIndels

    tmpFiles = {}
    for n in epitopeLens:
        tmpFasta = inFile.replace(".reformat.fasta",".tmp.%s.fasta"%(n))
        tmpFastaIndels = inFile.replace(".reformat.fasta",".tmp.%s.Indels.fasta"%(n))
        tmpFiles.update({n:tmpFasta})
        with open(tmpFasta, 'w') as outFile:
            for line in eps[n]:
                outFile.write(line)
        with open(tmpFastaIndels, 'w') as outFile:
            for line in epsIndels[n]:
                outFile.write(line)

    return(tmpFiles)

def ConstructAlleleHelper(s):
    return(s[:4].lower() + s[4:].capitalize())

def ConstructAlleles(hlas, FilePath, patID):
    '''
    Constructs the proper HLA input from HLA calls.

    :param hlas: list of HLA types for the Patient
    :return: list of normalized HLA identifiers for netMHCpan
    '''
    # TODO need a better way of verifying the format of the HLA alleles and matching in the list of those available...Some aren't working and should be...
    with open("%s/netMHCpanAlleles.txt"%(FilePath),'r') as alleles:
        allAlleles = [i.rstrip('\n').lower() for i in alleles.readlines()]
    
    hlas = [hla.lower() for hla in hlas if 'NA' not in hla]
    hlas = [i.replace("hla_","hla-") for i in hlas]
    hlas = [hla.replace("_","",1) for hla in hlas if 'NA' not in hla]
    hlas = [hla.replace("_",":",1) for hla in hlas if 'NA' not in hla]
    
    netMHCpanHLAS = []
    for hla in hlas:
        if len([h for h in allAlleles if hla == h]) == 1:
            netMHCpanHLAS.append(hla.replace("_", "").upper())
        elif len([h for h in allAlleles if hla.replace("_","")[0:-2] == h]) == 1:
            netMHCpanHLAS.append(hla.replace("_","")[0:-2].upper())
        elif len([h for h in allAlleles if hla.replace("_", "")[0:-4] == h]) == 1:
            netMHCpanHLAS.append(hla.replace("_", "")[0:-4].upper())
        elif len([h for h in allAlleles if hla.replace("_", "")[0:-5] == h]) == 1:
            netMHCpanHLAS.append(hla.replace("_", "")[0:-5].upper())
        else:
            sys.exit("ERROR: HLA type not found for %s %s" % (patID, hla))

    # for hla in hlas:
    #     if hla.replace("_","")[0:-2] in allAlleles:
    #         netMHCpanHLAS.append(hla.replace("_","")[0:-2].upper())
    #     elif hla.replace("_","")[0:-4] in allAlleles:
    #         netMHCpanHLAS.append(hla.replace("_", "")[0:-4].upper())
    #     elif hla.replace(" ", "") in allAlleles:
    #         netMHCpanHLAS.append(hla.replace("_", "").upper())
    #     else:
    #         sys.exit("ERROR: HLA type not found for %s %s" % (patID, hla))

    return(list(set(netMHCpanHLAS)))

def predict_neoantigens(FilePath, patName, inFile, hlasnormed, epitopeLens, netMHCpan, ELpred):
    '''
    Strips out all WILDTYPE and IMMEDIATE-STOPGAIN from fasta file.

    :param FilePath: Primary path of working directory
    :param patName: ID associated with a patient
    :param inFile: Fasta file with reformatted coding changes.
    :param hlas: HLA types for the patient.
    :param epitopeLens: List of epitope lengths to predict
    :param netMHCpan: Dictionary housing netMHCpan specific script locations and data. See README.md.
    :param ELpred: Logical for EL (true) or BA (false) predictions
    :return: netMHCpan predictions for each file.
    '''

    print("INFO: Predicting neoantigens for %s" % (patName))

    # Verify that the fasta file has information in it to avoid any errors thrown from netMHCpan
    checks = dict.fromkeys(inFile.keys())
    for n in inFile:
        cmd = "wc -l %s" % (inFile[n])
        pipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
        k = int(pipe.read().decode("utf-8").lstrip(" ").split(" ")[0])
        checks[n]=k

    epcalls = []
    for n in epitopeLens:
        if checks[n] > 0:
            output_file = 'tmp/%s.epitopes.%s.txt' % (patName, n)
            epcalls.append(output_file)
            with open(output_file, 'a') as epitope_pred:
                print("INFO: Running Epitope Predictions for %s on epitopes of length %s"%(patName,n))
                if ELpred:
                    cmd = [netMHCpan['netmhcpan'], '-l', str(n), '-a', ','.join(hlasnormed), '-f', inFile[n]]
                else:
                    cmd = [netMHCpan['netmhcpan'], '-BA', '-l', str(n), '-a', ','.join(hlasnormed), '-f', inFile[n]]
                netMHC_run = subprocess.Popen(cmd, stdout=epitope_pred, stderr=epitope_pred)
                netMHC_run.wait()
        else:
            print("INFO: Skipping Sample! No peptides to predict for %s" % (patName))

    print("INFO: Predictions complete for %s on epitopes of length %s" % (patName, n))

    return(epcalls)

def predict_neoantigensWT(FilePath, patName, inFile, hlasnormed, epitopeLens, netMHCpan):
    '''
    Strips out all WILDTYPE and IMMEDIATE-STOPGAIN from fasta file.

    :param FilePath: Primary path of working directory
    :param patName: ID associated with a patient
    :param inFile: Fasta file with reformatted coding changes.
    :param hlas: HLA types for the patient.
    :param epitopeLens: List of epitope lengths to predict
    :param netMHCpan: Dictionary housing netMHCpan specific script locations and data. See README.md.
    :return: netMHCpan predictions for each file.
    '''

    # Verify that the fasta file has information in it to avoid any errors thrown from netMHCpan
    checks = dict.fromkeys(inFile.keys())
    for n in inFile:
        cmd = "wc -l %s" % (inFile[n])
        pipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
        k = int(pipe.read().decode("utf-8").lstrip(" ").split(" ")[0])
        checks[n]=k

    epcalls = []
    for n in epitopeLens:
        if checks[n] > 0:
            output_file = '%s%s.wildtype.epitopes.%s.txt' % (FilePath, patName, n)
            epcalls.append(output_file)

            if os.path.isfile(output_file)==False:
                print("INFO: Predicting neoantigens for %s" % (patName))

                with open(output_file, 'a') as epitope_pred:
                    print("INFO: Running Epitope Predictions for %s on epitopes of length %s"%(patName,n))
                    cmd = [netMHCpan['netmhcpan'], '-BA', '-l', str(n), '-a', ','.join(hlasnormed), '-f', inFile[n]]
                    netMHC_run = subprocess.Popen(cmd, stdout=epitope_pred, stderr=epitope_pred)
                    netMHC_run.wait()
                    print("INFO: Predictions complete for %s on epitopes of length %s" % (patName, n))
            else:
                print("INFO: Neoantigen predictions already complete for %s epitopes of length %s" % (patName, n))
        else:
            print("INFO: Skipping Sample! No peptides to predict for %s" % (patName))

    return(epcalls)
