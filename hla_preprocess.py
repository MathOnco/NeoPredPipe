#!/usr/bin/env python

'''
@author: Eszter Lakatos & Ryan Schenck, (e.lakatos@qmul.ac.uk & ryan.schenck@univ.ox.ac.uk)
'''

import sys
import os
import glob


def readInHLAwinners(hladir):
    '''
    Transforms polysolver output to a single table line for hla file

    :param hladir: The directory which polysolver was outputted to
    :return: An array of HLA allele predictions, with NA for duplicates, as defined by NeoPredPipe
    '''
    hlaFileName = hladir.rstrip('/')+'/winners.hla.txt'
    with open(hlaFileName, 'r') as hlafile:
        hlaList = []
        for line in hlafile.readlines():
            lineArray = line.rstrip('\n').split('\t')
            hlaList.append(lineArray[1])
            if lineArray[2]!=lineArray[1]:
                hlaList.append(lineArray[2])
            else:
                hlaList.append('NA')
    return(hlaList)

def composeHLAFile(allHLAdir):
    outFileName = allHLAdir+'/hlatypes.txt'
    hlas = dict()
    with open(outFileName, 'w') as outFile:
        outFile.write('Patient\tHLA-A_1\tHLA-A_2\tHLA-B_1\tHLA-B_2\tHLA-C_1\tHLA-C_2\n')
        sampleList = glob.glob(allHLAdir+'/*')
        for sample in sampleList:
            hlaDir = sample
            sampleID = sample.split('/')[-1]
            hlaList = readInHLAwinners(hlaDir)
            outFile.write(sampleID+'\t'+ ('\t').join(hlaList)+'\n')
            hlas.update({sampleID: hlaList})
    return(hlas)


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
