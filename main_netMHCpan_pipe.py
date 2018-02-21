#!/usr/bin/env python

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
from vcf_manipulate import convert_to_annovar, annovar_annotation, get_coding_change,\
    predict_neoantigens, ReformatFasta, MakeTempFastas, ConstructAlleles
from postprocessing import DigestIndSample, AppendDigestedEps

def Parser():
    # get user variables
    parser = argparse.ArgumentParser()
    parser.add_argument("-E", "--epitopes", dest="epitopes", nargs='+', type=int, default=[8, 9, 10],
                        help="Epitope lengths for predictions. Default: 8 9 10")
    parser.add_argument("-l", dest="cleanLog", default=True, action='store_false', help="Specifies whether to delete the ANNOVAR log file. Default: True. Note: Use for debugging.")
    parser.add_argument("-d", dest="deleteInt", default=True, action='store_false', help="Specified whether to delete intermediate files created by program. Default: True. Note: Set flag to resume job.")
    parser.add_argument("-r", "--cleanrun", dest="makeitclean", default=False, action='store_true', help="Specify this alone with no other options to clean-up a run. Be careful that you mean to do this!!")
    parser.add_argument('-p', "--preponly", dest="preponly", default=False, action='store_true', help="Prep files only without running neoantigen predictions. The prediction step takes the most time.")
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
    postProcess.add_argument("-m", dest="checkPeptides", default=False, action='store_true',
                             help="Specifies whether to perform check if predicted epitopes match any normal peptide. If set to True, output is added as a column to neoantigens file. Requires PeptideMatch specified in usr_paths.ini. Default=False")
    postProcess.add_argument("-t", dest="buildSumTable", default=True, action='store_false', help="Flag to turn off a neoantigen burden summary table. Default=True.")

    Options = parser.parse_args()  # main user args
    if Options.makeitclean:
        CleanUp(Options)
        sys.exit("Process Complete")
    if not Options.vcfdir or not Options.hlafile or not Options.OutputDir:
        parser.error("Some of the required arguments were not provided. Please check required arguments.")

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

    def __init__(self, FilePath, patID, vcfFile, hla, annovar, netmhcpan, pepmatchPaths, Options):
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
        self.appendedEpitopes = None
        self.regionsPresent = None
        self.ProcessAnnovar(FilePath, annovar)
        self.hlasnormed = ConstructAlleles(self.hla, FilePath, self.patID)

        if Options.preponly:
            print("INFO: Input files prepared and completed for %s" % (self.patID))
        else:
            self.callNeoantigens(FilePath, netmhcpan, Options)
            if Options.postprocess:
                self.digestIndSample(pepmatchPaths, Options)

    def ProcessAnnovar(self, FilePath, annovar):
        # Prepare ANNOVAR input files
        if os.path.isfile("avready/"+self.patID+'.avinput'):
            print("INFO: ANNOVAR Ready files for %s already present."%(self.patID))
            self.avReadyFile = "avready/"+self.patID+'.avinput'
        else:
            self.avReadyFile = convert_to_annovar(FilePath, self.patID, self.vcfFile, annovar)

        # Prepare ANNOVAR annotated input files
        if os.path.isfile("avannotated/"+self.patID+'.avannotated.exonic_variant_function'):
            print("INFO: ANNOVAR Annotation files for %s already present." % (self.patID))
            self.annotationReady = "avannotated/"+self.patID+'.avannotated.exonic_variant_function'
        else:
            self.annotationReady = annovar_annotation(FilePath, self.patID, self.avReadyFile, annovar)

        # Get Coding Change
        if os.path.isfile("fastaFiles/"+self.patID+'.fasta'):
            print("INFO: Coding change fasta files for %s already present." % (self.patID))
            self.fastaChange = "fastaFiles/"+self.patID+'.fasta'
        else:
            self.fastaChange = get_coding_change(FilePath, self.patID, self.annotationReady, annovar)

    def callNeoantigens(self, FilePath, netmhcpan, Options):
        if os.path.isfile("fastaFiles/"+self.patID+'.reformat.fasta'):
            print("INFO: Coding change fasta files %s has already been reformatted." % (self.patID))
            self.fastaChangeFormat = "fastaFiles/"+self.patID+'.reformat.fasta'
        else:
            self.fastaChangeFormat = ReformatFasta(self.fastaChange)

        # Make tmp files.
        i = 0
        pepTmp = {}
        for n in Options.epitopes:
            if os.path.isfile("fastaFiles/%s.tmp.%s.fasta"%(self.patID,n)):
                pepTmp.update({n:"fastaFiles/%s.tmp.%s.fasta"%(self.patID,n)})
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
            if os.path.isfile("tmp/%s.epitopes.%s.txt" % (self.patID,n)):
                epTmp.append("tmp/%s.epitopes.%s.txt" % (self.patID,n))
                print("INFO: Epitope prediction files %s have already been created for netMHCpan length %s." % (self.patID,n))
                i += 1
                if i == len(Options.epitopes):
                    self.epcalls = epTmp
        if i!=len(Options.epitopes):
            self.epcalls = predict_neoantigens(FilePath, self.patID, self.peptideFastas, self.hlasnormed , Options.epitopes, netmhcpan)

    def digestIndSample(self, pmPaths, Options):
        if self.epcalls != []:
            checkPeptides=True
            self.digestedEpitopes = DigestIndSample(self.epcalls, self.patID, checkPeptides, pmPaths)
            self.appendedEpitopes, self.regionsPresent = AppendDigestedEps(self.digestedEpitopes, self.patID, self.annotationReady, self.avReadyFile, Options)

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

def FinalOut(sampleClasses, pepmatchPaths, Options):
    if Options.OutputDir[len(Options.OutputDir)-1]=="/":
        pass
    else:
        Options.OutputDir = Options.OutputDir + "/"

    if Options.includeall:
        outFile = Options.OutputDir + Options.outName + ".neoantigens.unfiltered.txt"
    else:
        outFile = Options.OutputDir + Options.outName + ".neoantigens.txt"

    outTable = Options.OutputDir + Options.outName + ".neoantigens.summarytable.txt"

    summaryTable = []
    with open(outFile, 'w') as pentultimateFile:
        if Options.includeall==True:
            for i in range(0, len(sampleClasses)):
                if sampleClasses[i].appendedEpitopes is not None:

                    pentultimateFile.write('\n'.join(sampleClasses[i].appendedEpitopes) + '\n')

                    for line in sampleClasses[i].appendedEpitopes:
                        if '<=' in line:
                            summaryTable.append(line)
        else:
            for i in range(0, len(sampleClasses)):
                if sampleClasses[i].appendedEpitopes is not None:
                    for line in sampleClasses[i].appendedEpitopes:
                        if '<=' in line:
                            pentultimateFile.write(line+"\n")
                            summaryTable.append(line)

    summaries = {}
    for z in range(0, len(sampleClasses)):
        if sampleClasses[z].appendedEpitopes is not None:
            total_count=0
            wbind=0
            sbind=0
            if Options.colRegions is not None:

                # Prep counts for multiregion data
                region_count=[0 for k in Options.colRegions*3]

                # Final Counters for each subtype of neoantigen
                clonal=0
                Wclonal = 0
                Sclonal = 0
                subclonal=0
                Wsubclonal = 0
                Ssubclonal=0
                shared=0
                Wshared = 0
                Sshared = 0

                regionsPesent = sampleClasses[z].regionsPresent
                overallRegions = Counter(regionsPesent)

            for line in sampleClasses[z].appendedEpitopes:
                # Counter to assess clonality of neoantigen
                r = 0
                rw = 0
                rs = 0

                if '<=' in line:
                    total_count+=1
                    if Options.colRegions is not None:
                        regions = [int(g) for g in line.split('\t')[1:len(Options.colRegions) + 1]]

                        for s in range(0,len(Options.colRegions)):
                            if regions[s] > 0:
                                region_count[s*3]+=1
                                r +=1

                    if "<=\tWB" in line:
                        wbind +=1

                        if Options.colRegions is not None:
                            regions = [int(g) for g in line.split('\t')[1:len(Options.colRegions) + 1]]
                            for s in range(0, len(Options.colRegions)):
                                if regions[s] > 0:
                                    region_count[(s * 3)+1] += 1
                                    rw +=1

                    elif "<=\tSB" in line:
                        sbind+=1

                        if Options.colRegions is not None:
                            regions = [int(g) for g in line.split('\t')[1:len(Options.colRegions) + 1]]
                            for s in range(0, len(Options.colRegions)):
                                if regions[s] > 0:
                                    region_count[(s * 3)+2] += 1
                                    rs+=1

                if Options.colRegions is not None:
                    if r == overallRegions['+']:
                        clonal += 1
                    elif r == 1:
                        subclonal += 1
                    elif r > 1 and overallRegions['+'] > 2:
                        shared += 1
                    elif r==0:
                        pass
                    else:
                        sys.exit('Error with mutliregion counter.')

                    if rw == overallRegions['+']:
                        Wclonal += 1
                    elif rw == 1:
                        Wsubclonal += 1
                    elif rw > 1 and overallRegions['+'] > 2:
                        Wshared += 1
                    elif rw==0:
                        pass
                    else:
                        sys.exit('Error with mutliregion counter.')

                    if rs == overallRegions['+']:
                        Sclonal += 1
                    elif rs == 1:
                        Ssubclonal += 1
                    elif rs > 1 and overallRegions['+'] > 2:
                        Sshared += 1
                    elif rs==0:
                        pass
                    else:
                        sys.exit('Error with mutliregion counter.')

            if Options.colRegions is not None:
                summaries.update({sampleClasses[z].patID:{'Total':total_count,'WB':wbind,'SB':sbind,
                                                        'Regions':region_count, 'Clonal':clonal, 'Subclonal':subclonal, 'Shared':shared,
                                                        'clonal_w':Wclonal, 'clonal_s':Sclonal, 'subclonal_w':Wsubclonal, 'subclonal_s':Ssubclonal,
                                                        'shared_w':Wshared,'shared_s':Sshared}})
            else:
                summaries.update({sampleClasses[z].patID:{'Total':total_count,'WB':wbind,'SB':sbind}})

    with open(outTable, 'w') as finalFile:
        if Options.colRegions is not None:
            header = ['Sample','Total','Total_WB','Total_SB','\t'.join(["Total_Region_%s"%(n) for n in range(0,len(Options.colRegions))]),
                      '\t'.join(["Total_WB_Region_%s" % (n) for n in range(0, len(Options.colRegions))]), '\t'.join(["Total_SB_Region_%s"%(n) for n in range(0,len(Options.colRegions))]),
                      'Clonal','Subclonal','Shared','Clonal_WB','Clonal_SB','Subclonal_WB','Subclonal_SB','Shared_WB','Shared_SB']

            finalFile.write('\t'.join(header) + '\n')

            for patient in summaries:
                line = [patient, summaries[patient]['Total'], summaries[patient]['WB'], summaries[patient]['SB'],
                        '\t'.join([str(summaries[patient]['Regions'][i]) for i in range(0,len(region_count), 3)]), '\t'.join([str(summaries[patient]['Regions'][i]) for i in range(1,len(region_count), 3)]),
                        '\t'.join([str(summaries[patient]['Regions'][i]) for i in range(2,len(region_count), 3)]),
                        summaries[patient]['Clonal'],summaries[patient]['Subclonal'],summaries[patient]['Shared'],
                        summaries[patient]['clonal_w'],summaries[patient]['clonal_s'],summaries[patient]['subclonal_w'],
                        summaries[patient]['subclonal_s'],summaries[patient]['shared_w'],summaries[patient]['shared_s']
                        ]
                line = [str(i) for i in line]
                finalFile.write('\t'.join(line) + '\n')
        else:
            header = ['Sample','Total','Total_WB','Total_SB']

            finalFile.write('\t'.join(header) + '\n')

            for patient in summaries:
                line = [patient, summaries[patient]['Total'], summaries[patient]['WB'], summaries[patient]['SB']]
                line = [str(i) for i in line]
                finalFile.write('\t'.join(line) + '\n')

    print("INFO: Summary Tables Complete.")

def CleanUp(Options):
    if Options.cleanLog or Options.makeitclean:
        try:
            os.remove('logforannovarNeoPredPipe.txt')
        except OSError:
            pass
        try:
            os.remove('logForPeptideMatch.tmp')
        except OSError:
            pass
    if (Options.deleteInt or Options.makeitclean) and not Options.preponly:
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
    localpath = os.path.abspath(__file__).replace('main_netMHCpan_pipe.py', '')  # path to scripts working directory
    Config = configparser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    annPaths = ConfigSectionMap(Config.sections()[0], Config)  # get annovar script paths
    netMHCpanPaths = ConfigSectionMap(Config.sections()[1], Config)  # get annovar script paths
    try:
        pepmatchPaths = ConfigSectionMap(Config.sections()[2], Config)  # get PeptideMatch paths
    except IndexError as e:
        pepmatchPaths = None

    Options = Parser()
    print("INFO: Begin.")

    allFiles, hlas = PrepClasses(localpath, Options)

    # Prepare samples
    t = []
    for patFile in allFiles:
        patname = patFile.replace('.vcf', '').split("/")[len(patFile.replace('.vcf', '').split("/")) - 1]
        t.append(Sample(localpath, patname, patFile, hlas[patname], annPaths, netMHCpanPaths, pepmatchPaths, Options))

    if Options.preponly:
        print("INFO: Complete.")
	print("INFO: Preprocessed intermediary files are in avready, avannotated and fastaFiles. If you wish to perform epitope prediction, run the pipeline again without the --preponly flag, intermediary files will be automatically detected.")
    else:
        if Options.postprocess:
            FinalOut(t, pepmatchPaths, Options)
            print("INFO: Complete")
        else:
            print("INFO: Complete")

    CleanUp(Options)

if __name__=="__main__":
    main()
