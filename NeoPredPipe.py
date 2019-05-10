#!/usr/bin/env python

'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
Contributions from: Eszter Lakatos
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
from vcf_manipulate import convert_to_annovar, annovar_annotation, get_coding_change,\
    ReformatFasta, MakeTempFastas
from predict_binding import predict_neoantigens
from hla_preprocess import ConstructAlleles, ConstructAlleles_typeII, composeHLAFile, composeHLA2File, readInHLAwinners
from postprocessing import DigestIndSample, AppendDigestedEps
from process_expression import GetExpressionFiles

def Parser():
    # get user variables
    parser = argparse.ArgumentParser()
    parser.add_argument("-E", "--epitopes", dest="epitopes", nargs='+', type=int, default=[8, 9, 10],
                        help="Epitope lengths for predictions. Default: 8 9 10")
    parser.add_argument("-l", dest="cleanLog", default=True, action='store_false', help="Specifies whether to delete the ANNOVAR log file. Default: True. Note: Use for debugging.")
    parser.add_argument("-d", dest="deleteInt", default=True, action='store_false', help="Specifies whether to delete intermediate files created by program. Default: True. Note: Set flag to resume job.")
    parser.add_argument("-r", "--cleanrun", dest="makeitclean", default=False, action='store_true', help="Specify this alone with no other options to clean-up a run. Be careful that you mean to do this!!")
    parser.add_argument('-p', "--preponly", dest="preponly", default=False, action='store_true',help="Prep files only without running neoantigen predictions. The prediction step takes the most time.")
    parser.add_argument("--EL", dest="ELpred", default=False, action='store_true',
        help="Flag to perform netMHCpan predictions with Eluted Ligand option (without the -BA flag). Please note that the output will NOT be compatible with downstream Recognition Potential analysis. Default=False (BA predictions)")
    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument("-I", dest="vcfdir", default=None, type=str,
                               help="Input vcf file directory location. Example: -I ./Example/input_vcfs/")
    requiredNamed.add_argument("-H", dest="hlafile", default=None, type=str, help="HLA file for vcf patient samples OR directory with patient-specific directories from running POLYSOLVER (see Readme).")
    requiredNamed.add_argument("-o", dest="OutputDir", default=None, type=str, help="Output Directory Path")
    requiredNamed.add_argument("-n", dest="outName", default="AllSamples", type=str, help="Name of the output file for neoantigen predictions")
    postProcess = parser.add_argument_group('Post Processing Options')
    postProcess.add_argument("-pp", dest="postprocess", default=True, action='store_false', help="Flag to perform post processing. Default=True.")
    postProcess.add_argument("-c", dest="colRegions", default=None, nargs='+',
                        help="Columns of regions within vcf that are not normal within a multiregion vcf file after the format field. Example: 0 is normal in test samples, tumor are the other columns. Program can handle different number of regions per vcf file.")
    postProcess.add_argument("-a", dest="includeall", default=False, action='store_true', help="Flag to not filter neoantigen predictions and keep all regardless of prediction value.")
    postProcess.add_argument("-m", dest="checkPeptides", default=False, action='store_true',
                             help="Specifies whether to perform check if predicted epitopes match any normal peptide. If set to True, output is added as a column to neoantigens file. Requires PeptideMatch specified in usr_paths.ini. Default=False")
    postProcess.add_argument("-x", "--expression", dest="expression", default=None, type=str,
                        help="RNAseq expression quantification file(s), if specified, expression information is added to output tables.")
    postProcess.add_argument("--expmulti", dest="expMultiregion", default=False, action='store_true',
                        help="Flag to specify if expression file(s) has information on multiple regions in multiple columns. Default=False.")
    postProcess.add_argument("-t", dest="buildSumTable", default=True, action='store_false', help="Flag to turn off a neoantigen burden summary table. Default=True.")

    Options = parser.parse_args()  # main user args
    Options.typeII = False
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
        self.ProcessAnnovar(FilePath, annovar, Options)
        if Options.typeII:
            self.hlasnormed = ConstructAlleles_typeII(self.hla, FilePath, self.patID)
        else:
            self.hlasnormed = ConstructAlleles(self.hla, FilePath, self.patID)

        if Options.preponly:
            print("INFO: Input files prepared and completed for %s" % (self.patID))
        else:
            self.callNeoantigens(FilePath, netmhcpan, Options)
            if Options.postprocess:
                self.digestIndSample(FilePath, pepmatchPaths, Options)

    def ProcessAnnovar(self, FilePath, annovar, Options):
        # Prepare ANNOVAR input files
        if os.path.isfile(Options.OutputDir+"avready/"+self.patID+'.avinput'):
            print("INFO: ANNOVAR Ready files for %s already present."%(self.patID))
            self.avReadyFile = Options.OutputDir+"avready/"+self.patID+'.avinput'
        else:
            self.avReadyFile = convert_to_annovar(Options.OutputDir, self.patID, self.vcfFile, annovar)

        # Prepare ANNOVAR annotated input files
        if os.path.isfile(Options.OutputDir+"avannotated/"+self.patID+'.avannotated.exonic_variant_function'):
            print("INFO: ANNOVAR Annotation files for %s already present." % (self.patID))
            self.annotationReady = Options.OutputDir+"avannotated/"+self.patID+'.avannotated.exonic_variant_function'
        else:
            self.annotationReady = annovar_annotation(Options.OutputDir, self.patID, self.avReadyFile, annovar)

        # Get Coding Change
        if os.path.isfile(Options.OutputDir+"fastaFiles/"+self.patID+'.fasta'):
            print("INFO: Coding change fasta files for %s already present." % (self.patID))
            self.fastaChange =  Options.OutputDir+"fastaFiles/"+self.patID+'.fasta'
        else:
            self.fastaChange = get_coding_change(Options.OutputDir, self.patID, self.annotationReady, annovar)

        if os.path.isfile(Options.OutputDir+"fastaFiles/"+self.patID+'.reformat.fasta'):
            print("INFO: Coding change fasta files %s has already been reformatted." % (self.patID))
            self.fastaChangeFormat = Options.OutputDir+"fastaFiles/"+self.patID+'.reformat.fasta'
        else:
            self.fastaChangeFormat = ReformatFasta(self.fastaChange)

    def callNeoantigens(self, FilePath, netmhcpan, Options):
        # Make tmp files.
        i = 0
        pepTmp = {}
        for n in Options.epitopes:
            if os.path.isfile(Options.OutputDir+"fastaFiles/%s.tmp.%s.fasta"%(self.patID,n)) and os.path.isfile("fastaFiles/%s.tmp.%s.fasta"%(self.patID,str(n)+'.Indels')):
                pepTmp.update({n:Options.OutputDir+"fastaFiles/%s.tmp.%s.fasta"%(self.patID,n)})
                pepTmp.update({str(n)+'.Indels':Options.OutputDir+"fastaFiles/%s.tmp.%s.fasta"%(self.patID,str(n)+'.Indels')})
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
            if os.path.isfile(Options.OutputDir+"tmp/%s.epitopes.%s.txt" % (self.patID,n)):
                epTmp.append(Options.OutputDir+"tmp/%s.epitopes.%s.txt" % (self.patID,n))
                i += 1
                print("INFO: Epitope prediction files %s have already been created for netMHCpan length %s." % (self.patID,n))
            if os.path.isfile(Options.OutputDir+"tmp/%s.epitopes.%s.txt" % (self.patID,str(n)+'.Indels')):
                epTmp.append(Options.OutputDir+"tmp/%s.epitopes.%s.txt" % (self.patID,str(n)+'.Indels'))
                i += 1
                print("INFO: Epitope prediction files %s have already been created for netMHCpan length %s." % (self.patID,str(n)+'.Indels'))
            if i == 2*len(Options.epitopes):
                self.epcalls = epTmp
        if i!=2*len(Options.epitopes):
            if i>0:
                os.system("rm "+Options.OutputDir+"tmp/"+self.patID+".epitopes.*.txt") # if doing predictions, remove existing files to ensure double predicting happens
            self.epcalls = predict_neoantigens(Options.OutputDir, self.patID, self.peptideFastas, self.hlasnormed , Options.epitopes, netmhcpan, Options.ELpred)

    def digestIndSample(self, FilePath, pmPaths, Options):
        if self.epcalls != []:
            toDigestSNVs = filter(lambda y: 'Indels.txt' not in y, self.epcalls)
            toDigestIndels = filter(lambda y: 'Indels.txt' in y, self.epcalls)
            if toDigestSNVs != []:
                self.digestedEpitopes = DigestIndSample(toDigestSNVs, self.patID, Options, pmPaths)
                self.appendedEpitopes, self.regionsPresent = AppendDigestedEps(FilePath, self.digestedEpitopes, self.patID, self.annotationReady, self.avReadyFile, Options)
            else:
                self.appendedEpitopes = None
                self.regionsPresent = None
            if toDigestIndels != []:
                self.digestedEpitopesIndels = DigestIndSample(toDigestIndels, self.patID, Options, pmPaths, True)
                self.appendedEpitopesIndels, self.regionsPresentIndels = AppendDigestedEps(FilePath, self.digestedEpitopesIndels, self.patID, self.annotationReady, self.avReadyFile, Options)
            else:
                self.appendedEpitopesIndels = None
                self.regionsPresentIndels = None

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
    if os.path.isfile(Options.hlafile):
        with open(Options.hlafile, 'r') as hlaFile:
            lines = hlaFile.readlines()
        for line in lines:
            line = line.rstrip('\n').split("\t")
            pat = line[0]
            del line[0]
            hlas.update({pat: line})
    elif os.path.isdir(Options.hlafile):
        if Options.typeII:
            hlas = composeHLA2File(Options.hlafile)
        else:
            hlas = composeHLAFile(Options.hlafile)
    else:
        sys.exit("Unable to locate HLA file or directory.")

    if not os.path.exists(Options.OutputDir):
        os.mkdir(Options.OutputDir)
    
    try:
        os.mkdir(Options.OutputDir+'avready')
        
    except OSError as e:
        print("INFO: Proper directory already exists. Continue.")

    try:
        os.mkdir(Options.OutputDir+'avannotated')
    except OSError as e:
        print("INFO: Proper directory already exists. Continue.")

    try:
        os.mkdir(Options.OutputDir+'fastaFiles')
    except OSError as e:
        print("INFO: Proper directory already exists. Continue.")

    try:
        os.mkdir(Options.OutputDir+'tmp')
    except OSError as e:
        print("INFO: Proper directory already exists. Continue.")

    return(allFiles, hlas)

def FinalOut(sampleClasses, Options, indelProcess=False):
    if indelProcess:
        epitopesToProcess = 'appendedEpitopesIndels' #Process this special set of predicted epitopes
        regionsToProcess = 'regionsPresentIndels'
        filePostFix = '.neoantigens.Indels'
    else:
        epitopesToProcess = 'appendedEpitopes'
        regionsToProcess = 'regionsPresent'
        filePostFix = '.neoantigens'

    if Options.includeall:
        outFile = Options.OutputDir + Options.outName + filePostFix + ".unfiltered.txt"
    else:
        outFile = Options.OutputDir + Options.outName + filePostFix + ".txt"

    outTable = Options.OutputDir + Options.outName + filePostFix + ".summarytable.txt"

    summaryTable = []
    with open(outFile, 'w') as pentultimateFile:
        if Options.includeall==True:
            for i in range(0, len(sampleClasses)):
                appendedEps = getattr(sampleClasses[i],epitopesToProcess)
                if appendedEps is not None:

                    pentultimateFile.write('\n'.join(appendedEps) + '\n')

                    for line in appendedEps:
                        if '<=' in line:
                            summaryTable.append(line)
        else:
            for i in range(0, len(sampleClasses)):
                appendedEps = getattr(sampleClasses[i],epitopesToProcess)
                if appendedEps is not None:
                    for line in appendedEps:
                        if '<=' in line:
                            pentultimateFile.write(line+"\n")
                            summaryTable.append(line)

    summaries = {}
    for z in range(0, len(sampleClasses)):
        appendedEps = getattr(sampleClasses[z],epitopesToProcess)
        if appendedEps is not None:
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

                regionsPesent = getattr(sampleClasses[z],regionsToProcess)
                overallRegions = Counter(regionsPesent)

            for line in appendedEps:
                # Counter to assess clonality of neoantigen
                r = 0
                rw = 0
                rs = 0

                if '<=' in line and not (Options.checkPeptides and line[-1]=='0'):
                    total_count+=1
                    if Options.colRegions is not None:
                        regions = [int(g) for g in line.split('\t')[1:len(Options.colRegions) + 1]]

                        for s in range(0,len(Options.colRegions)):
                            if regions[s] > 0:
                                region_count[s*3]+=1
                                r +=1

                    if "<=\tWB" in line or "<=WB" in line:
                        wbind +=1

                        if Options.colRegions is not None:
                            regions = [int(g) for g in line.split('\t')[1:len(Options.colRegions) + 1]]
                            for s in range(0, len(Options.colRegions)):
                                if regions[s] > 0:
                                    region_count[(s * 3)+1] += 1
                                    rw +=1

                    elif "<=\tSB" in line or "<=SB" in line:
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
            os.remove(Options.OutputDir+'logforannovarNeoPredPipe.txt')
        except OSError:
            pass
        try:
            os.remove(Options.OutputDir+'logForPeptideMatch.tmp')
        except OSError:
            pass
    if (Options.deleteInt or Options.makeitclean) and not Options.preponly:
        try:
            shutil.rmtree(Options.OutputDir+"avready/")
        except OSError as e:
            print("ERROR: Unable to clean intermediary files.")
            print(e)
        try:
            shutil.rmtree(Options.OutputDir+"avannotated/")
        except OSError as e:
            print("ERROR: Unable to clean intermediary files.")
            print(e)
        try:
            shutil.rmtree(Options.OutputDir+"fastaFiles/")
        except OSError as e:
            print("ERROR: Unable to clean intermediary files.")
            print(e)
        try:
            shutil.rmtree(Options.OutputDir+"tmp/")
        except OSError as e:
            print("ERROR: Unable to clean tmp files.")
            print(e)

def main():
    # Pull information about usr system files
    localpath = os.path.abspath(__file__).replace('NeoPredPipe.py', '')  # path to scripts working directory
    Options = Parser()
    Config = configparser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    annPaths = ConfigSectionMap(Config.sections()[0], Config)  # get annovar script paths
    annPaths['build'] = annPaths['gene_table'].split('/')[-1].split('_')[0] # get build version from name of gene table (hg19/hg38_refGene...)
    annPaths['gene_model'] = annPaths['gene_table'].split('/')[-1].split('_')[1].replace(".txt","") # get gene model from X_somethingGene.txt
    if annPaths['build'] not in ['hg19', 'hg38', 'hg18']:
        print("WARNING: Unexpected genome build detected in annovar reference files: %s. Please check path for 'gene_table'. Build hg19 is used for analysis."%(annPaths['build']))
        annPaths['build'] = 'hg19'
    else:
        print("INFO: Annovar reference files of build %s were given, using this build for all analysis."%(annPaths['build']))
    netMHCpanPaths = ConfigSectionMap(Config.sections()[1], Config)  # get netmhcpan paths
    if netMHCpanPaths['netmhcpan'].rstrip('\n').split('/')[-1]=='netMHCIIpan': #set typeII if netMHCIIpan is supplied instead of netMHCpan
        Options.typeII = True
    try:
        pepmatchPaths = ConfigSectionMap(Config.sections()[2], Config)  # get PeptideMatch paths
    except IndexError as e:
        pepmatchPaths = None

    #Make sure outputdir ends with '/'
    if Options.OutputDir[len(Options.OutputDir)-1]=="/":
        pass
    else:
        Options.OutputDir = Options.OutputDir + "/"

    #Check if PeptideMatch paths are provided, ignore -m if not
    if Options.checkPeptides and pepmatchPaths is None:
        print("WARNING: You chose to perform peptide match checking for epitopes, but did not provide paths for PeptideMatch. The pipeline will ignore the -m flag")
        Options.checkPeptides = False

    print("INFO: Begin.")

    allFiles, hlas = PrepClasses(localpath, Options)
    # Eliminate header line and tailing empty rows
    hlas.pop('', None)
    hlas.pop('Patient', None)

    # Expression related options
    if Options.expression is not None:
        Options.allExpFiles = GetExpressionFiles(Options)

    # Check VCF and HLA
    assert len(allFiles) > 0, "No input vcf files detected. Perhaps they are compressed?"
    if len(allFiles)>len(hlas.keys()):
        print("WARNING: Less samples are detected in HLA file than in VCF folder. Only samples included in HLA file will be processed.")

    # Prepare samples
    t = []
    for patFile in allFiles:
        patname = patFile.replace('.vcf', '').split("/")[len(patFile.replace('.vcf', '').split("/")) - 1]
        if patname in hlas.keys():
            t.append(Sample(localpath, patname, patFile, hlas[patname], annPaths, netMHCpanPaths, pepmatchPaths, Options))

    if len(t)<len(hlas.keys()):
    	print("WARNING: Not all samples in HLA file have matching VCF files. Please check that HLA file is tab-separated and sample names match exactly with .vcf file names. Only matching samples will be included in analysis and output tables.")

    if Options.preponly:
        print("INFO: Complete.")
        print("INFO: Preprocessed intermediary files are in avready, avannotated and fastaFiles. If you wish to perform epitope prediction, run the pipeline again without the --preponly flag, intermediary files will be automatically detected.")
    else:
        if Options.postprocess:
            FinalOut(t, Options)
            FinalOut(t, Options, True)
            print("INFO: Complete")
        else:
            print("INFO: Complete")

    CleanUp(Options)

if __name__=="__main__":
    main()
