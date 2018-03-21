#!/usr/bin/env python

'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
'''

import sys
import os
import glob
import argparse
try:
    import ConfigParser as configparser # for python 2
except:
    import configparser # for python 3
import shutil
from collections import Counter
from NeoClass import Neoantigen
from NeoAlign import Aligner

def Parser():
    # get user variables
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dirty", dest="Dirty", default=True, action='store_false', help="Flag to keep intermediate files. Default: False. Note: Use for debugging.")
    parser.add_argument("-o", '--neoreco_out', dest="neorecoOut", default="", type=str, help="Output Directory.")
    parser.add_argument("-a", '--midpoint', dest='a', default=1., type=float, help="Midpoint parameter of the logistic function, alignment score threshold.")
    parser.add_argument("-k", '--slope', dest='k', default=1., type=float, help="Slope parameter of the logistic function")
    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument("-i", '--neopred_in', dest="neoPredIn", default=None, type=str,help="Input neoantigen predictions. Example: -I ./AllSamples.neoantigens.txt")
    requiredNamed.add_argument("-f", '--fastas', dest="fastaDir", default="fastaFiles/", type=str,help="Fasta files directory associated with the input.")
    Options = parser.parse_args()  # main user args
    if not Options.neoPredIn:
        parser.error("Some of the required arguments were not provided. Please check required arguments.")
    if len(Options.neorecoOut)!=0:
        if Options.neorecoOut[len(Options.neorecoOut)-1]!='/':
            Options.neorecoOut = Options.neorecoOut+'/'
    if len(Options.fastaDir)!=0:
        if Options.fastaDir[len(Options.fastaDir)-1]!='/':
            Options.fastaDir = Options.fastaDir+'/'
    return(Options)

def ConfigSectionMap(section, Config):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                print("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

class StandardPreds:
    '''
    Holds standard predicitons from the other pipeline
    '''
    def __init__(self, Options):
        self.filename = Options.neoPredIn
        self.fastaPath = Options.fastaDir
        self.samples = []
        self.hlas = None
        self.fastas = None
        self.preds = None

    def load(self):
        '''
        Loads the data class of neoantigen predictions to get the right information.

        :return:
        '''
        with open(self.filename, 'r') as inFile:
            lines = [line.replace('\n','') for line in inFile.readlines()]
        self.preds = lines
        self.samples = list(set([line.split('\t')[0] for line in lines]))
        self.hlas = {sam:[] for sam in self.samples}
        self.fastas = {sam:'%s%s.reformat.fasta'%(self.fastaPath,sam) for sam in self.samples}

        for line in lines:
            sam = line.split('\t')[0]
            hla = line.split('\t')[11]
            if hla not in self.hlas[sam]:
                self.hlas[sam].append(hla)

def GetWildTypePredictions():
    pass


def main():
    # Pull information about usr system files
    localpath = os.path.abspath(__file__).replace('NeoRecoPo.py', '')  # path to scripts working directory
    Config = configparser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    Options = Parser()

    tmpOut = '%s/%s/NeoRecoTMP/'%(Options.neorecoOut, localpath)
    if os.path.isdir(tmpOut)==False:
        os.system('mkdir %s'%(tmpOut))

    preds = StandardPreds(Options) # Create instance of StandardPreds Class
    preds.load() # Load the neoantigen predictions data

    if Options.Dirty:
        os.system('rm -r %s'%(tmpOut))


if __name__=="__main__":
    main()