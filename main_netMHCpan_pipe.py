#!/usr/bin/env python

import sys
import os
import glob
import argparse
import ConfigParser
import shutil
from vcf_manipulate import convert_to_annovar, annovar_annotation, get_coding_change

def Parser():
    # get user variables
    parser = argparse.ArgumentParser()
    parser.add_argument("-E", "--epitopes", dest="epitopes", nargs='+', type=int, default=[8, 9, 10],
                        help="Epitope lengths for predictions. Default: 8 9 10")
    parser.add_argument("-m", dest="multiregion", default=False, action='store_true',
                        help="Specifies if the vcf is a multiregion sample. Default: False.")
    parser.add_argument("-c", dest="colRegions", default=None, nargs="+",
                        help="Columns of regions within vcf that are not normal multiregion vcf file. 0 is normal in test samples. Can handle different number of regions per vcf file.")
    parser.add_argument("-l", dest="cleanLog", default=True, action='store_false', help="Specifies whether to delete the ANNOVAR log file. Default: True. Note: Use for debugging.")
    parser.add_argument("-d", dest="deleteInt", default=True, action='store_false', help="Specified whether to delete intermediate files created by program. Default: True. Note: Set flag to resume job.")
    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument("-I", dest="vcfdir", default=None, type=str,
                               help="Input vcf file directory location. Example: -I ./Example/input_vcfs/")
    requiredNamed.add_argument("-H", dest="hlafile", default=None, type=str, help="HLA file for vcf patient samples.")
    requiredNamed.add_argument("-o", dest="OutputDir", default=None, type=str, help="Output Directory Path")

    # parser.add_option('-r', dest='permute', default=False, action='store_true', help='Permute sequences [Default: %default]')
    Options = parser.parse_args()  # main user args
    if not Options.vcfdir or not Options.hlafile or not Options.OutputDir:
        parser.error("Some of the required arguments were not provided. Please check required arguments.")
    if Options.multiregion and Options.colRegions is None:
        parser.error("-m requires -c to be specified.")

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

class Sample():

    def __init__(self, FilePath, patID, vcfFile, hla, annovar, netmhcpan):
        self.patID = patID
        self.vcfFile = vcfFile
        self.hla = hla
        self.avReadyFile=None
        self.annotationReady=None
        self.ProcessAnnovar(FilePath, annovar)

    def loadvcfs(self):
        with open(self.vcfFile, 'r') as inputFile:
            muts = [line.rstrip("\n") for line in inputFile.readlines()]

    def ProcessAnnovar(self, FilePath, annovar):
        # Prepare ANNOVAR input files
        if os.path.isfile(FilePath+"avready/"+self.patID+'.avinput'):
            print("INFO: ANNOVAR Ready files for %s already present."%(self.patID))
            self.avReadyFile = FilePath+"avready/"+self.patID+'.avinput'
        else:
            self.avReadyFile = convert_to_annovar(FilePath, self.patID, self.vcfFile, annovar)

        # Prepare ANNOVAR annotated input files
        if os.path.isfile(FilePath+"avannotated/"+self.patID+'.avannotated.exonic_variant_function'):
            print("INFO: ANNOVAR Annotation files for %s already present." % (self.patID))
            self.annotationReady = FilePath+"avannotated/"+self.patID+'.avannotated.exonic_variant_function'
        else:
            self.annotationReady = annovar_annotation(FilePath, self.patID, self.avReadyFile, annovar)

        # Get Coding Change
        if os.path.isfile(FilePath+"fastaFiles/"+self.patID+'.fasta'):
            print("INFO: Coding change fasta files for %s already present." % (self.patID))
            self.fastaChange = FilePath+"fastaFiles/"+self.patID+'.fasta'
        else:
            self.fastaChange = get_coding_change(FilePath, self.patID, self.annotationReady, annovar)


def PrepClasses(FilePath, Options):
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

    try:
        os.mkdir('avready')
    except OSError as e:
        print("INFO: Proper directory already exist. Continue.")

    try:
        os.mkdir('avannotated')
    except OSError as e:
        print("INFO: Proper directory already exist. Continue.")

    try:
        os.mkdir('fastaFiles')
    except OSError as e:
        print("INFO: Proper directory already exist. Continue.")

    return(allFiles, hlas)

def CleanUp(Options):
    if Options.cleanLog:
        try:
            os.remove('logforannovarNeoPredPipe.txt')
        except OSError:
            pass
    if Options.deleteInt:
        try:
            shutil.rmtree("avready/")
        except OSError as e:
            print("ERROR: Unable to clean intermediary files.")
            print(e)
        try:
            shutil.rmtree("avannotated/")
        except OSError as e:
            print("ERROR: Unable to clean intermediary files.")
            print(e)
        try:
            shutil.rmtree("fastaFiles/")
        except OSError as e:
            print("ERROR: Unable to clean intermediary files.")
            print(e)

def main():
    print("Info: Begin.")
    # Pull information about usr system files
    localpath = os.path.abspath(__file__).rstrip('main_netMHCpan_pipe.py')  # path to scripts working directory
    Config = ConfigParser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    annPaths = ConfigSectionMap(Config.sections()[0], Config)  # get annovar script paths
    netMHCpanPaths = ConfigSectionMap(Config.sections()[1], Config)  # get annovar script paths
    Options = Parser()

    allFiles, hlas = PrepClasses(localpath, Options)

    # Prepare samples
    t = []
    for patFile in allFiles:
        patname = patFile.rstrip(".vcf").split("/")[len(patFile.rstrip(".vcf").split("/")) - 1]
        t.append(Sample(localpath, patname, patFile, hlas[patname], annPaths, netMHCpanPaths))

    CleanUp(Options)


if __name__=="__main__":
    main()
