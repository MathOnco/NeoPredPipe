#!/usr/bin/env python

'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
'''
import sys
import os
import argparse
from collections import OrderedDict
try:
    import ConfigParser as configparser # for python 2
except:
    import configparser # for python 3

def Parser():
    # get user variables
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument("-o", '--output', dest="output", default="", type=str, help="Output Directory.")
    requiredNamed.add_argument("-i", '--neopred_in', dest="neoPredIn", default=None, type=str,help="Input neoantigen predictions, must be filtered on binding affinity. Example: -i ./AllSamples.neoantigens.txt")
    requiredNamed.add_argument("-w", '--wt_neo_table', dest="wt_neo_table", default=None, type=str,help="Intermediary file created by NeoRecoPo, this contains the WT and neoantigen binding affinities.")
    requiredNamed.add_argument("-r", '--RecoPreds', dest="recopreds", default=None, type=str,help="Final output of NeoRecoPo containing the Recognition Potentials.")
    requiredNamed.add_argument("-s", '--SampleTypes', dest="SamTypes", default=None)
    Options = parser.parse_args()  # main user args

    if not Options.neoPredIn or not Options.wt_neo_table:
        parser.error("Some of the required arguments were not provided. Please check required arguments.")

    return(Options)

def loadPreds(fileName):
    refDict = {}
    with open(fileName, 'r') as inputFile:
        lines = [line.replace('\n', '').split('\t') for line in inputFile.readlines()]

    for i, line in enumerate(lines):
        if len(line) < 26:
            line.append('')
            line.append('')
            lines[i] = line

    if lines[0][0].startswith("Sample"):
        header = lines[0][0]
    else:
        regionNums = 0
        for i in range(1, len(lines[0])):
            if lines[0][i].startswith('line'):
                break
            else:
                regionNums += 1
        header = '\t'.join(
            ['Sample', '\t'.join(['R' + str(i) for i in range(1, regionNums + 1)]), 'Line', 'Chr', 'AllelePos', 'Ref',
             'Alt', 'Symbol', 'Pos', 'hla', 'peptide', 'core', 'Of', 'Gp', 'Gl', 'lp', 'll', 'Icore', 'Identity',
             'Score', 'BA', 'Rank', 'Candidate', 'BindLevel']).split('\t')
    data = [dict(zip(header, v)) for v in lines[0:]]
    for neo in data:
        id = '.'.join([neo['Sample'], neo['Identity'], str(len(neo['peptide']))])
        refDict.update({id: neo})

    return(header, refDict)

def loadRecos(fileName):
    refDict = {}
    with open(fileName, 'r') as inputFile:
        lines = [line.replace('\n', '').split('\t') for line in inputFile.readlines()]

    header = lines[0]

    data = [dict(zip(header, v)) for v in lines[0:]]
    for neo in data:
        id = '.'.join([neo['Sample'], neo['Mutation'], str(len(neo['MutantPeptide']))])
        refDict.update({id: neo})

    return (header, refDict)

def loadWTMutTable(fileName):
    refDict = {}
    with open(fileName, 'r') as inputFile:
        lines = [line.replace('\n', '').split('\t') for line in inputFile.readlines()]

    header = lines[0]

    data = [dict(zip(header, v)) for v in lines[0:]]
    for neo in data:
        id = '.'.join([neo['Sample'], neo['MUTATION_ID'], str(len(neo['MT.PEPTIDE']))])
        refDict.update({id: neo})

    return (header, refDict)

def BuildMasterTable(localpath, Options):
    predHeader, preds = loadPreds(Options.neoPredIn)
    recoHeader, reco = loadRecos(Options.recopreds)
    wtHeader, wtTab = loadWTMutTable(Options.wt_neo_table)

    for neo in preds:
        a = preds[neo]
        b = reco[neo]
        c = wtTab[neo]



def main():
    print("INFO: Begin.")
    # Pull information about usr system files
    localpath = os.path.abspath(__file__).replace('NeoRecoPoBlend.py', '')  # path to scripts working directory

    Config = configparser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    Options = Parser()

    BuildMasterTable(localpath, Options)


if __name__=="__main__":
    main()