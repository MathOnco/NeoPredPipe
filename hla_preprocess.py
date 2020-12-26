#!/usr/bin/env python

'''
@author: Eszter Lakatos & Ryan Schenck, (e.lakatos@qmul.ac.uk & ryan.schenck@univ.ox.ac.uk)
'''

import sys
import os
import glob
import itertools
import re



def readInHLAwinners(hladir):
    '''
    Transforms polysolver output to a single table line for hla file

    :param hladir: The directory which polysolver was outputted to
    :return: An array of HLA allele predictions, with NA for duplicates, as defined by NeoPredPipe
    '''
    hlaFileName = hladir.rstrip('/')+'/winners.hla.txt'
    if not os.path.isfile(hlaFileName):
    	return(None)
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
            if hlaList != None:
            	outFile.write(sampleID+'\t'+ ('\t').join(hlaList)+'\n')
            	hlas.update({sampleID: hlaList})
    return(hlas)

def processHLAminerFile(hlaFileName):
    '''
    Transforms HLAminer output (which is really not a csv) into a list of best hits of HLA2 haplotypes (<=2 preds per allele)

    :param hladir: The directory which HLAminer was outputted to
    :return: A filtered file with the filtered alleles
    '''
    pattern = r'D(\w+)\*\d+:\d+[PNQG]?[,:]\d+[.:]\d+,\d+\.\d+e?[-,]\d+[,.]\d+e?-?(\d+)?[.,]\d+\.?(\d+)?'

    filteredfile = open(hlaFileName.replace('HPRA.csv','processed.txt'),'w')
    hlaMainList = []
    with open(hlaFileName, 'r') as hlafile: 
        lines = hlafile.readlines()
        for line in lines:
            if re.search(pattern, line):
                hla = line.strip('\t').split(',')[0]
                hlaMain = hla.split(':')[0]
                if not hlaMain in hlaMainList:
                    hlaMainList.append(hlaMain)
                    filteredfile.write(hla+'\n')
    filteredfile.close()

def readInHLA2hlaminer(hlaFileName):
    '''
    Transforms pre-processed HLAminer output to a single table line for hla file

    :param hladir: The directory which HLAminer was outputted to
    :return: An array of HLA allele predictions
    '''
    hlaProcessedFileName = hlaFileName.replace('HPRA.csv','processed.txt')
    if not os.path.isfile(hlaProcessedFileName):
        processHLAminerFile(hlaFileName)
    with open(hlaProcessedFileName, 'r') as hlafile:
        hlaList = []
        for line in hlafile.readlines():
            hla = line.rstrip('\n')
            #filter out DRA and DRB2-DRB9, change this line if it is not required
            if not re.search(r'DR[AB][^1]',hla):
                hlaList.append(hla)
    return(hlaList)

def readInHLA2hlahd(hlaFileName):
    '''
    Transforms pre-processed HLA-HD output to a single table line for hla file

    :param hladir: The directory which HLA-HD was outputted to
    :return: An array of HLA allele predictions
    '''
    with open(hlaFileName, 'r') as hlafile:
        hlaList = []
        for line in hlafile.readlines():
            lineArray = line.rstrip('\n').split('\t')
            # This line should be modified to include A-C or DBR2-9 or other predictions
            if re.search(r'D[PQR][AB]1', lineArray[0]):
                for hla in lineArray[1:]:
                    if hla!='-' and hla!='Not typed':
                        hlaList.append(hla[4:])
    return(hlaList)

def composeHLA2File(allHLAdir):
    outFileName = allHLAdir+'/hlatypes.txt'
    hlas = dict()
    with open(outFileName, 'w') as outFile:
        sampleList = glob.glob(allHLAdir+'/*')
        for sample in sampleList:
            hladir = sample
            sampleID = sample.split('/')[-1]
            if os.path.isfile(hladir.rstrip('/')+'/result/'+sampleID+'_final.result.txt'):
                hlaList = readInHLA2hlahd(hladir.rstrip('/')+'/result/'+sampleID+'_final.result.txt')
            elif os.path.isfile(hladir.rstrip('/')+'/HLAminer_HPRA.csv'):
                hlaList = readInHLA2hlaminer(hladir.rstrip('/')+'/HLAminer_HPRA.csv')
            else:
                hlaList = None
            if hlaList is not None:
                outFile.write(sampleID+'\t'+ ('\t').join(hlaList)+'\n')
                hlas.update({sampleID: hlaList})
    return(hlas)


def ConstructAlleleHelper(s):
    return(s[:4].lower() + s[4:].capitalize())

def ConstructAllelesOld(hlas, FilePath, patID):
    '''
    Constructs the proper HLA input from HLA calls.

    :param hlas: list of HLA types for the Patient
    :return: list of normalized HLA identifiers for netMHCpan
    '''
    # Outdated version only kept for historical reasons for now!
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

    return(list(set(netMHCpanHLAS)))

def ConstructAlleles(hlas, FilePath, patID):
    '''
    Constructs the proper HLA input from HLA calls.

    :param hlas: list of HLA types for the Patient
    :return: list of normalized HLA identifiers for netMHCpan
    '''

    with open("%s/netMHCpanAlleles.txt"%(FilePath),'r') as alleles:
        allAlleles = [i.rstrip('\n').upper() for i in alleles.readlines()]
    
    hlas = [hla.upper() for hla in hlas]
    hlas = [hla for hla in hlas if 'NA' not in hla]
    if 'H-2' not in hlas[0] and 'H2' not in hlas[0]:
        hlasWithSuffix = [hla for hla in hlas if hla[-1] in 'NQSALC']
        hlas = [hla.rstrip('NQSALC') for hla in hlas]
    hlas = [i.replace("HLA_","HLA-") for i in hlas]
    hlas = ['_'.join(hla.split('_')[:3]) for hla in hlas] #Strip 5th-8th digits as they do not modify protein structure
    hlas = [hla.replace("_","",1) for hla in hlas]
    hlas = [hla.replace("_",":",1) for hla in hlas]

    if len(hlasWithSuffix)>0:
        print("WARNING: Expression status indicating suffix found for %s : %s.\n Keep in mind when analysing results!" % (patID, ', '.join(hlasWithSuffix)))
    
    netMHCpanHLAS = []
    for hla in hlas:
        if hla in allAlleles:
            netMHCpanHLAS.append(hla)
        else:
            sys.exit("ERROR: HLA type not found for %s %s" % (patID, hla))

    return(list(set(netMHCpanHLAS)))

def ConstructAlleles_typeII(hlas, FilePath, patID):
    '''
    Constructs the proper HLA input from HLA calls.

    :param hlas: list of HLA types for the Patient
    :return: list of normalized HLA identifiers for netMHCpan
    '''
    with open("%s/netMHCpanAlleles2.txt"%(FilePath),'r') as alleles:
        allAlleles = [i.rstrip('\n').upper() for i in alleles.readlines()]
    
    hlas = [hla.upper() for hla in hlas]
    hlas = [hla.rstrip("PNQG") for hla in hlas if 'NA' not in hla] #potential characters at the end that have to be omitted
    hlas = [''.join(hla.split(':')[:2]) for hla in hlas]

    hlas_drb = []; hlas_dpa = []; hlas_dpb = []; hlas_dqa = []; hlas_dqb = []
    for hla in hlas:
        if 'DRB' in hla:
            hlas_drb.append(hla.replace("*","_",1))
        elif 'DPA' in hla:
            hlas_dpa.append(hla.replace("*","",1))
        elif 'DPB' in hla:
            hlas_dpb.append(hla.replace("*","",1))
        elif 'DQA' in hla:
            hlas_dqa.append(hla.replace("*","",1))
        elif 'DQB' in hla:
            hlas_dqb.append(hla.replace("*","",1))

    hlas_dp = ['HLA-'+'-'.join(hla) for hla in list(itertools.product(hlas_dpa, hlas_dpb))]
    hlas_dq = ['HLA-'+'-'.join(hla) for hla in list(itertools.product(hlas_dqa, hlas_dqb))]

    netMHCpanHLAS = []
    for hla in hlas_drb:
        if hla in allAlleles:
            netMHCpanHLAS.append(hla)
        else:
            print("ERROR (non-fatal): HLA type not found for %s : %s. Continuing with other hlas." % (patID, hla))
    for hla in hlas_dp + hlas_dq:
        if hla in allAlleles:
            netMHCpanHLAS.append(hla)
        else:
            print("WARNING: HLA type combination not found for %s: %s. It will be omitted in the analysis." % (patID, hla))
    

    return(list(set(netMHCpanHLAS)))
