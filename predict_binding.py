#!/usr/bin/env python

'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
Contributions from: Eszter Lakatos
'''

import sys
import os
import subprocess

def predict_neoantigens(FilePath, patName, inFile, hlasnormed, epitopeLens, netMHCpan, Options):
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
    for n in inFile:
        if checks[n] > 0:
            output_file = FilePath+'tmp/%s.epitopes.%s.txt' % (patName, n)
            epcalls.append(output_file)
            with open(output_file, 'a') as epitope_pred:
                print("INFO: Running Epitope Predictions for %s on epitopes of length %s"%(patName,n))
                if Options.ELpred:
                    cmd = [netMHCpan['netmhcpan'], '-l', str(n).split('.')[0], '-a', ','.join(hlasnormed), '-f', inFile[n]]
                elif Options.typeII:
                    cmd = [netMHCpan['netmhcpan'], '-length', str(n).split('.')[0], '-a', ','.join(hlasnormed), '-f', inFile[n]]
                else:
                    cmd = [netMHCpan['netmhcpan'], '-BA', '-l', str(n).split('.')[0], '-a', ','.join(hlasnormed), '-f', inFile[n]]
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
