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
def convert_to_annovar(FilePath, patName, inFile, annovar, manualproc=False):
    '''
    :param FilePath: Primary path of working directory
    :param patName: ID associated with a patient
    :param inFile: VCF file for neoantigens to be predicted from
    :param annovar: Dictionary housing annovar specific script locations. See README.md.
    :param manualproc: Boolean flag indicating if manual processing should be done instead of annovar's convert2annovar
    :return: Name of the ANNOVAR ready vcf file for annovar_annnotation()
    '''
    print("INFO: Running convert2annovar.py on %s" % (inFile))

    # specify name for output file
    outDir = FilePath+"avready/"
    annovar_out_ready = outDir + patName + '.avinput'

    if manualproc:
        # process vcf file manually for avinput format
        with open(inFile, 'r') as vcfFile:
            vcfLines = vcfFile.readlines()
        with open(annovar_out_ready, 'w') as processedFile:
            for line in vcfLines:
                if line[0]=='#':
                    processedFile.write(line.rstrip()+'\n')
                else:
                    linespl = line.split('\t')
                    vcf_file_chr = linespl[0]
                    vcf_file_start = linespl[1]
                    vcf_file_end = linespl[1]
                    vcf_file_ref = linespl[3]
                    vcf_file_alt = linespl[4]
                    vcf_file_qual = linespl[5]
                    if( len( vcf_file_ref ) > len( vcf_file_alt ) ): # fix to work with deletions (tested on Sterlka2 vcf output)
                        deletion_length = len( vcf_file_ref ) - len( vcf_file_alt )
                        vcf_file_end = str( int( vcf_file_start ) + deletion_length )
                    processedFile.write('\t'.join([ vcf_file_chr, vcf_file_start, vcf_file_end, vcf_file_ref, vcf_file_alt, '.', vcf_file_qual, '.', line.rstrip() ]) + '\n') #frequency and depth info is filled with placeholder
        print('INFO: Manual VCF Conversion Process complete %s'%(inFile))

    else:
        with open(FilePath+"logforannovarNeoPredPipe.txt", 'a') as logFile:
            cmd = ['perl', annovar['convert2annovar'], '-format', 'vcf4', inFile, '-outfile', annovar_out_ready, '-allsample',
                    '-includeinfo', '-withfreq', '-comment']
            runconvert = subprocess.Popen(cmd, stdout=logFile, stderr=logFile)
            runconvert.wait()
        print('INFO: ANNOVAR VCF Conversion Process complete %s'%(inFile))

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
    outDir = FilePath+"avannotated/"
    annovar_out_ready = outDir + patName + '.avannotated'

    with open(FilePath+"logforannovarNeoPredPipe.txt", 'a') as logFile:
        cmd = ['perl', annovar['annotatevariation'], '-out', annovar_out_ready, '-build', annovar['build'], '-dbtype', annovar['gene_model'], inFile, annovar['humandb'], '--comment']
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
    outDir = FilePath+"fastaFiles/"
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
            if 'wildtype' not in seq_record.id.lower() and 'immediate-stopgain' not in seq_record.id.lower() and 'silent' not in seq_record.id.lower() and 'startloss' not in seq_record.id.lower():
                # TODO: Add a regex expression to extract the position since it's variable with versions of ANNOVAR
                # TODO: Regex code for this might be r"\w*((?i)position;\d+;(?-i))\W*"
                try:
                    pos = int(seq_record.id.replace(";;",";").split(";")[5].split('-')[0])-1
                except ValueError:
                    try:
                        pos = int(seq_record.id.replace(";;",";").split(";")[6].split('-')[0])-1
                    except IndexError:
                        sys.exit("ERROR: Could not process fasta line in reformatted fasta: %s" % (seq_record.id))

                if 'dup' in seq_record.id.lower() or 'del' in seq_record.id.lower() or 'ins' in seq_record.id.lower() or 'from;*;to;' in seq_record.id.lower() or 'fs' in seq_record.id: #if mutation is potentially frameshift
                    miniseq = ExtractSeq(seq_record, pos, n, True) # flag for frameshift
                    mySeqsIndels.append(">"+seq_record.id[0:100]+"\n"+miniseq+"\n")
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
        tmpFiles.update({str(n)+'.Indels':tmpFastaIndels})
        with open(tmpFasta, 'w') as outFile:
            for line in eps[n]:
                outFile.write(line)
        with open(tmpFastaIndels, 'w') as outFile:
            for line in epsIndels[n]:
                outFile.write(line)

    return(tmpFiles)
