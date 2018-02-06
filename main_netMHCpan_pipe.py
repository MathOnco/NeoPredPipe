#!/usr/bin/env python

import sys
import os
import glob
import argparse
import ConfigParser
import shutil
from vcf_manipulate import convert_to_annovar, annovar_annotation, get_coding_change,\
    predict_neoantigens, ReformatFasta, MakeTempFastas
from postprocessing import DigestIndSample, AppendDigestedEps

def Parser():
    # get user variables
    parser = argparse.ArgumentParser()
    parser.add_argument("-E", "--epitopes", dest="epitopes", nargs='+', type=int, default=[8, 9, 10],
                        help="Epitope lengths for predictions. Default: 8 9 10")
    parser.add_argument("-l", dest="cleanLog", default=True, action='store_false', help="Specifies whether to delete the ANNOVAR log file. Default: True. Note: Use for debugging.")
    parser.add_argument("-d", dest="deleteInt", default=True, action='store_false', help="Specified whether to delete intermediate files created by program. Default: True. Note: Set flag to resume job.")
    parser.add_argument("-r", "--cleanrun", dest="makeitclean", default=False, action='store_true', help="Specify this alone with no other options to clean-up a run. Be careful that you mean to do this!!")
    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument("-I", dest="vcfdir", default=None, type=str,
                               help="Input vcf file directory location. Example: -I ./Example/input_vcfs/")
    requiredNamed.add_argument("-H", dest="hlafile", default=None, type=str, help="HLA file for vcf patient samples.")
    requiredNamed.add_argument("-o", dest="OutputDir", default=None, type=str, help="Output Directory Path")
    requiredNamed.add_argument("-n", dest="outName", default="AllSamples", type=str, help="Name of the output file for neoantigen predictions")
    postProcess = parser.add_argument_group('Post Processing Options')
    postProcess.add_argument("-pp", dest="postprocess", default=True, action='store_false', help="Flag to perform post processing. Default=True.")
    postProcess.add_argument("-c", dest="colRegions", default=None, nargs="+",
                        help="Columns of regions within vcf that are not normal within a multiregion vcf file after the format field. Example: 0 is normal in test samples, tumor are the other columns. Program can handle different number of regions per vcf file.")
    postProcess.add_argument("-a", dest="includeall", default=False, action='store_true', help="Flag to not filter neoantigen predictions and keep all regardless of prediction value.")
    postProcess.add_argument("-t", dest="buildSumTable", default=True, action='store_false', help="Flag to turn off summary table.")

    Options = parser.parse_args()  # main user args
    if Options.makeitclean:
        CleanUp(Options)
        sys.exit("Process Complete")
    if not Options.vcfdir or not Options.hlafile or not Options.OutputDir:
        parser.error("Some of the required arguments were not provided. Please check required arguments.")
    if Options.postprocess:
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
    '''
    Use this to run and execute the pipeline on an individual patient vcf file.
    '''

    def __init__(self, FilePath, patID, vcfFile, hla, annovar, netmhcpan, Options):
        self.patID = patID
        self.vcfFile = vcfFile
        self.hla = hla
        self.avReadyFile=None
        self.annotationReady=None
        self.fastaChange=None
        self.fastaChangeFormat=None
        self.peptideFastas = None # Will be a dictionary of tmp files for predictions
        self.epcalls = None
        self.digestedEpitopes = None
        self.ProcessAnnovar(FilePath, annovar)
        self.callNeoantigens(FilePath, netmhcpan, Options)
        if Options.postprocess:
            self.digestIndSample(Options)

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

    def callNeoantigens(self, FilePath, netmhcpan, Options):
        if os.path.isfile(FilePath+"fastaFiles/"+self.patID+'.reformat.fasta'):
            print("INFO: Coding change fasta files %s has already been reformatted." % (self.patID))
            self.fastaChangeFormat = FilePath+"fastaFiles/"+self.patID+'.reformat.fasta'
        else:
            self.fastaChangeFormat = ReformatFasta(self.fastaChange)

        # Make tmp files.
        i = 0
        pepTmp = {}
        for n in Options.epitopes:
            if os.path.isfile(FilePath+"fastaFiles/%s.tmp.%s.fasta"%(self.patID,n)):
                pepTmp.update({n:FilePath+"fastaFiles/%s.tmp.%s.fasta"%(self.patID,n)})
                print("INFO: Tmp fasta files %s has already been created for netMHCpan length %s." % (self.patID,n))
                i+=1
                if i == len(Options.epitopes):
                    self.peptideFastas = pepTmp
        if i != len(Options.epitopes):
            self.peptideFastas = MakeTempFastas(self.fastaChangeFormat, Options.epitopes)

        # Predict neoantigens
        i = 0
        epTmp = []
        for n in Options.epitopes:
            if os.path.isfile(FilePath + "tmp/%s.epitopes.%s.txt" % (self.patID,n)):
                epTmp.append(FilePath + "tmp/%s.epitopes.%s.txt" % (self.patID,n))
                print("INFO: Epitope prediction files %s have already been created for netMHCpan length %s." % (self.patID,n))
                i += 1
                if i == len(Options.epitopes):
                    self.epcalls = epTmp
        if i!=len(Options.epitopes):
            self.epcalls = predict_neoantigens(FilePath, self.patID, self.peptideFastas, self.hla, Options.epitopes, netmhcpan)

    def digestIndSample(self, Options):
        self.digestedEpitopes = DigestIndSample(self.epcalls, self.patID)
        self.appendedEpitopes = AppendDigestedEps(self.digestedEpitopes, self.patID, self.annotationReady, self.avReadyFile, Options)

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
        print("INFO: Proper directory already exists. Continue.")

    try:
        os.mkdir('avannotated')
    except OSError as e:
        print("INFO: Proper directory already exists. Continue.")

    try:
        os.mkdir('fastaFiles')
    except OSError as e:
        print("INFO: Proper directory already exists. Continue.")

    try:
        os.mkdir('tmp')
    except OSError as e:
        print("INFO: Proper directory already exists. Continue.")

    return(allFiles, hlas)

def CleanUp(Options):
    if Options.cleanLog or Options.makeitclean:
        try:
            os.remove('logforannovarNeoPredPipe.txt')
        except OSError:
            pass
    if Options.deleteInt or Options.makeitclean:
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
        try:
            shutil.rmtree("tmp/")
        except OSError as e:
            print("ERROR: Unable to clean tmp files.")
            print(e)

def main():
    # Pull information about usr system files
    localpath = os.path.abspath(__file__).rstrip('main_netMHCpan_pipe.py')  # path to scripts working directory
    Config = ConfigParser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    annPaths = ConfigSectionMap(Config.sections()[0], Config)  # get annovar script paths
    netMHCpanPaths = ConfigSectionMap(Config.sections()[1], Config)  # get annovar script paths
    Options = Parser()
    print("INFO: Begin.")

    allFiles, hlas = PrepClasses(localpath, Options)

    # Prepare samples
    t = []
    for patFile in allFiles:
        patname = patFile.rstrip(".vcf").split("/")[len(patFile.rstrip(".vcf").split("/")) - 1]
        t.append(Sample(localpath, patname, patFile, hlas[patname], annPaths, netMHCpanPaths, Options))

    CleanUp(Options)

if __name__=="__main__":
    main()
