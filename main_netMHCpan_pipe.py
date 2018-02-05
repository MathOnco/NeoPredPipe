#!/usr/bin/env python

import sys
import os
import glob
import argparse
import ConfigParser

def Parser():
    # get user variables
    parser = argparse.ArgumentParser()
    parser.add_argument("-E", "--epitopes", dest="epitopes", nargs='+', type=int, default=[8, 9, 10],
                        help="Epitope lengths for predictions. Default: 8 9 10")

    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument("-I", dest="vcfdir", default=None, type=str,
                               help="Input vcf file directory location. Example: -I ./Example/input_vcfs/")
    requiredNamed.add_argument("-H", dest="hlafile", default=None, type=str, help="HLA file for vcf patient samples.")
    requiredNamed.add_argument("-o", dest="OutputDir", default=None, type=str, help="Output Directory Path")
    # parser.add_option('-r', dest='permute', default=False, action='store_true', help='Permute sequences [Default: %default]')
    Options = parser.parse_args()  # main user args

    if not Options.vcfdir or not Options.hlafile or not Options.OutputDir:
        sys.exit("Some of the required arguments were not provided. Please check required arguments.")

    return(Options)

def ConfigSectionMap(section):
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

class Sample():

    def __init__(self, patID, vcfFile, hla):
        self.patID = patID
        self.vcfFile = vcfFile
        self.hla = hla
        self.muts = self.loadvcfs()

    def loadvcfs(self):
        with open(self.vcfFile, 'r') as inputFile:
            muts = [line.rstrip("\n") for line in inputFile.readlines()]

def PrepClasses(Options):
    if Options.vcfdir[len(Options.vcfdir)-1] != "/":
        Options.vcfdir+="/"
    getFiles = Options.vcfdir + "*.vcf"
    allFiles = glob.glob(getFiles)
    if all([os.path.isfile(f) for f in allFiles]):
        pass
    else:
        getFiles = Options.vcfdir + "/*.vcf"
        allFiles = glob.glob(getFiles)
        if all([os.path.isfile(f) for f in allFiles]):
            pass
        else:
            sys.exit("Unable to locate vcf files.")

    hlas = dict()
    try:
        with open(Options.hlafile, 'r') as hlaFile:
            lines = hlaFile.readlines()
    except:
        sys.exit("Unable to locate HLA file.")
    for line in lines:
        line = line.rstrip('\n').split("\t")
        pat = line[0]
        del line[0]
        hlas.update({pat: line})

    t = []
    for patFile in allFiles:
        patname = patFile.rstrip(".vcf").split("/")[len(patFile.rstrip(".vcf").split("/"))-1]
        t.append(Sample(patname, patFile, hlas[patname]))
    return(t)

if __name__=="__main__":
    # Pull information about usr system files
    localpath = os.path.abspath(__file__).rstrip('main_netMHCpan_pipe.py')  # path to scripts working directory
    Config = ConfigParser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    annPaths =ConfigSectionMap(Config.sections()[0]) # get annovar script paths
    Options = Parser()

    samples = PrepClasses(Options)

