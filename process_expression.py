#!/usr/bin/env python

'''
@author: Eszter Lakatos, e.lakatos@qmul.ac.uk & Ryan Schenck, ryan.schenck@univ.ox.ac.uk
'''

import os
import glob

def GetExpressionFiles(Options):
    '''
    Fetches expression file(s) from user input

    :param Options.expression: either a single expression file, or a folder containing expression files in the format <sampleX>.tsv
    :return: filenames in
    '''
    if os.path.exists(Options.expression):
        if os.path.isdir(Options.expression):
            getFiles = Options.vcfdir + "*.vcf"
            allExpFiles = glob.glob(Options.expression + "/*.tsv")
        else:
            allExpFiles = Options.expression
    else:
        print("WARNING: Unable to locate expression file(s), please check file path. Flag -x is ignored.")
        allExpFiles = None
        Options.expression = None

    return(allExpFiles)

def BuildGeneIDTable(FilePath,idType):
    '''
    Links RefSeq IDs (coming from Annovar) to Ensembl or USCS IDs

    :param idType: 'ensembl' (gene or transcript) or 'uscs'
    :return: dictionary to map NM IDs to other IDs
    '''
    idDict = {'ensembl_gene':0, 'ensembl_transcript':1, 'uscs':4, 'name':3, 'refseq':2}
    idInd = idDict[idType]

    with open(FilePath+'mart_table_hg38_unique.txt', 'r') as martTable:
        lines = martTable.readlines()
        header = lines[0]
        lines = lines[1:]
        idTable = {line.split('\t')[2]:line.split('\t')[idInd].rstrip('\n') for line in lines}

    return(idTable)

def BuildExpTable(expFile, multiregion):
    '''
    Reads in expression data into dictionary

    :param expFile: expression file path
    :return: dictionary of gene IDs and expression values
    '''
    with open(expFile, 'r') as expression:
        lines = expression.readlines()
        if not lines[0].split('\t')[1][0].isdigit(): #first line is header and needs to be excluded
            lines = lines[1:]

        testID = lines[0].split('\t')[0]
        idType = 'refseq'
        if testID[0:4]=='ENST':
            idType = 'ensembl_transcript'
        elif testID[0:4]=='ENSG':
            idType = 'ensembl_gene'
        elif testID[0:2]=='uc':
            idType = 'uscs'

        if multiregion:
            expTable = {line.split('\t')[0].split('.')[0]:','.join(line.rstrip('\n').split('\t')[1:])  for line in lines}
        else:
            expTable = {line.split('\t')[0].split('.')[0]:line.split('\t')[1].rstrip('\n') for line in lines}

    return(idType, expTable)



