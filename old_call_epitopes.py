#!/usr/bin/env python2.7

# Send the fasta files to the various predictors

import sys
from Bio import SeqIO
import glob
import subprocess
import os

# input fasta file directory
arg1 = sys.argv[1]

# output directory
arg2 = sys.argv[2]

# length of epitope
arg3 = int(sys.argv[3])


# get fasta files into a dictionary
def get_files_for_prediction(arg1):
    # get the full path to the fasta files
    arg1_full = os.path.abspath(arg1) + '/'

    # get all of the fasta files in the folder
    for_prediction = glob.glob(arg1_full + '*.fasta')

    # get file and path and it's name into a dictionary
    # key = filename
    # value = path + file
    for_epitope_pred = dict()
    for item in for_prediction:
        file_name = item[item.rfind('/') + 1:item.rfind('.')]
        file_name = file_name.split('_new_head')
        file_name = file_name[0]
        for_epitope_pred.update({file_name: item})

    return for_epitope_pred


# creates a temporary fasta file with non-wildtype sequences only per vcf file
def non_wildtype_seq(for_epitope_pred, allele_length):
    '''STRIPS OUT:
    WILD-TYPE
    IMMEDIATE-STOPGAIN'''

    tmp_fastas = dict()
    '''
    for item in for_epitope_pred:
        fasta_file = 'tmp/' + item + '.fasta'

        tmp_fastas.update({item:fasta_file})

        with open(fasta_file, 'w') as temp_fasta:
            for seq_record in SeqIO.parse(for_epitope_pred[item], 'fasta'):
                if 'wildtype' in seq_record.id.lower():
                    pass
                elif 'immediate-stopgain' in seq_record.id.lower():
                    pass
                else:
                    temp_fasta.write('>'+seq_record.id + '\n')
                    temp_fasta.write(str(sub_seq) + '\n')
    '''

    if allele_length == 10:
        for item in for_epitope_pred:
            fasta_file = 'tmp/' + item + '.fasta'

            tmp_fastas.update({item: fasta_file})

            with open(fasta_file, 'w') as temp_fasta:
                for seq_record in SeqIO.parse(for_epitope_pred[item], 'fasta'):
                    if 'wildtype' in seq_record.id.lower():
                        pass
                    elif 'immediate-stopgain' in seq_record.id.lower():
                        pass
                    else:
                        mut_pos = str(seq_record.id).split('(')
                        mut_pos = mut_pos[1].split(';')
                        mut_pos = int(mut_pos[1])

                        if (len(seq_record) - mut_pos) >= 9 and mut_pos > 9:
                            # Captures all those that are flanked by proper numbers of residues
                            sub_seq = seq_record.seq[mut_pos - 10:mut_pos + 9]
                            temp_fasta.write('>' + seq_record.id + '\n')
                            temp_fasta.write(str(sub_seq) + '\n')

                        elif mut_pos <= 9:
                            # captures all those at the beginning of the sequence
                            sub_seq = seq_record.seq[0:mut_pos + 9]
                            temp_fasta.write('>' + seq_record.id + '\n')
                            temp_fasta.write(str(sub_seq) + '\n')

                        elif (len(seq_record) - mut_pos) < 9:
                            # captures the end of the sequence
                            sub_seq = seq_record.seq[mut_pos - 10:len(seq_record)]
                            temp_fasta.write('>' + seq_record.id + '\n')
                            temp_fasta.write(str(sub_seq) + '\n')

                        else:
                            pass

    elif allele_length == 9:
        for item in for_epitope_pred:
            fasta_file = 'tmp/' + item + '.fasta'

            tmp_fastas.update({item: fasta_file})

            with open(fasta_file, 'w') as temp_fasta:
                for seq_record in SeqIO.parse(for_epitope_pred[item], 'fasta'):
                    if 'wildtype' in seq_record.id.lower():
                        pass
                    elif 'immediate-stopgain' in seq_record.id.lower():
                        pass
                    else:
                        mut_pos = str(seq_record.id).split('(')
                        mut_pos = mut_pos[1].split(';')
                        mut_pos = int(mut_pos[1])

                        if (len(seq_record) - mut_pos) >= 8 and mut_pos > 8:
                            # Captures all those that are flanked by proper numbers of residues
                            sub_seq = seq_record.seq[mut_pos - 9:mut_pos + 8]
                            temp_fasta.write('>' + seq_record.id + '\n')
                            temp_fasta.write(str(sub_seq) + '\n')

                        elif mut_pos <= 8:
                            # captures all those at the beginning of the sequence
                            sub_seq = seq_record.seq[0:mut_pos + 8]
                            temp_fasta.write('>' + seq_record.id + '\n')
                            temp_fasta.write(str(sub_seq) + '\n')

                        elif (len(seq_record) - mut_pos) < 8:
                            # captures the end of the sequence
                            sub_seq = seq_record.seq[mut_pos - 9:len(seq_record)]
                            temp_fasta.write('>' + seq_record.id + '\n')
                            temp_fasta.write(str(sub_seq) + '\n')

                        else:
                            pass

    elif allele_length == 8:
        for item in for_epitope_pred:
            fasta_file = 'tmp/' + item + '.fasta'

            tmp_fastas.update({item: fasta_file})

            with open(fasta_file, 'w') as temp_fasta:
                for seq_record in SeqIO.parse(for_epitope_pred[item], 'fasta'):
                    if 'wildtype' in seq_record.id.lower():
                        pass
                    elif 'immediate-stopgain' in seq_record.id.lower():
                        pass
                    else:
                        mut_pos = str(seq_record.id).split('(')
                        mut_pos = mut_pos[1].split(';')
                        mut_pos = int(mut_pos[1])

                        if (len(seq_record) - mut_pos) >= 7 and mut_pos > 7:
                            # Captures all those that are flanked by proper numbers of residues
                            sub_seq = seq_record.seq[mut_pos - 8:mut_pos + 7]
                            print(len(sub_seq))
                            temp_fasta.write('>' + seq_record.id + '\n')
                            temp_fasta.write(str(sub_seq) + '\n')

                        elif mut_pos <= 7:
                            # captures all those at the beginning of the sequence
                            sub_seq = seq_record.seq[0:mut_pos + 7]
                            # print len(sub_seq)
                            temp_fasta.write('>' + seq_record.id + '\n')
                            temp_fasta.write(str(sub_seq) + '\n')

                        elif (len(seq_record) - mut_pos) < 7:
                            # captures the end of the sequence
                            sub_seq = seq_record.seq[mut_pos - 8:len(seq_record)]
                            temp_fasta.write('>' + seq_record.id + '\n')
                            temp_fasta.write(str(sub_seq) + '\n')

                        else:
                            pass

    return tmp_fastas


# netMHC predictor
def netMHC(tmp_fastas, arg2, allele_length):
    # allele_file = '/Users/Ryan/Desktop/TEST/alleles/netMHC_alleles/netMHC-4.0_HLA_supertype_representative.txt'
    allele_file = '/Users/Ryan/Desktop/TEST/alleles/netMHC_alleles/netMHC_HLA.txt'

    lengths = ['8', '9', '10']
    output_directory = arg2
    for item in tmp_fastas:
        output_file = output_directory + item + '_netMHC_predictions.txt'
        print 'Running predictions for %s' % (tmp_fastas[item])
        with open(output_file, 'a') as epitope_pred:
            with open(allele_file, 'r') as alleles:
                for allele in alleles:
                    print 'Running predictions for %s' % (allele.rstrip('\n'))
                    print 'Epitope length: %s' % (allele_length)
                    allele = allele.rstrip('\n')
                    # run program for each allele
                    netMHC_run = subprocess.Popen(
                        ['netMHC', '-l', str(allele_length), '-a', allele, '-f', tmp_fastas[item]], stdout=epitope_pred)
                    netMHC_run.wait()

# clean up temp files
def cleanup(tmp_fastas):
    # delete all temp fastas
    for tmp_file in tmp_fastas:
        os.remove(tmp_fastas[tmp_file])

def main(arg1, arg2, arg3):
    os.makedirs('tmp')

    # get fasta files to execute predictions in dictionary
    for_epitope_pred = get_files_for_prediction(arg1)

    # create custom temporary fasta files for epitope predictors
    tmp_fastas = non_wildtype_seq(for_epitope_pred, arg3)

    # perform netMHC predictions
    netMHC(tmp_fastas, arg2, arg3)

    # clean up temp environment
    cleanup(tmp_fastas)
    # remove temporary directory
    os.rmdir('tmp')


main(arg1, arg2, arg3)