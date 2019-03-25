#!/usr/bin/env python

'''
@author: Eszter Lakatos, e.lakatos@qmul.ac.uk & Ryan Schenck, ryan.schenck@univ.ox.ac.uk
'''

import sys
import os

def BuildGeneIDTable(idType):
    '''
    Links RefSeq IDs (coming from Annovar) to Ensembl or USCS IDs

    :param idType: 'ensembl' (gene or transcript) or 'uscs'
    :return: dictionary to map NM IDs to other IDs
    '''
    idDict = {'ensembl_gene':0, 'ensembl_transcript':1, 'uscs':4, 'name':3, 'refseq':2}
    idInd = idDict[idType]

    with open('mart_table_hg38_unique.txt', 'r') as martTable:
        lines = martTable.readlines()
        header = lines[0]
        lines = lines[1:]
        idTable = {line.split('\t')[2]:line.split('\t')[idInd].rstrip('\n') for line in lines}

    return(idTable)

def BuildExpTable(expFile):
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

        expTable = {line.split('\t')[0].split('.')[0]:float(line.split('\t')[1].rstrip('\n')) for line in lines}

    return(idType, expTable)



