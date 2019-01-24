#!/usr/bin/env python

'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
'''

import sys
import os
import argparse
try:
    import ConfigParser as configparser # for python 2
except:
    import configparser # for python 3
from Bio import SeqIO
import pickle
from NeoClass import Neoantigen
from NeoAlign import Aligner

from StandardPredsClass import StandardPreds

def Parser():
    # get user variables
    parser = argparse.ArgumentParser()
    parser.add_argument("-id", "--indel", dest="Indels", default=False, action='store_true', help="Treat inputted antigen predictions as indels/frameshifts. Wild-type peptide sequences will not be used for RecognitionPotential calculation.")
    parser.add_argument("-d", "--dirty", dest="Dirty", default=True, action='store_false', help="Flag to keep intermediate files. Default: False. Note: Use for debugging.")
    parser.add_argument("-o", '--neoreco_out', dest="neorecoOut", default="", type=str, help="Output Directory.")
    parser.add_argument("-a", '--midpoint', dest='a', default=26.0, type=float, help="Midpoint parameter of the logistic function, alignment score threshold.")
    parser.add_argument("-k", '--slope', dest='k', default=4.86936, type=float, help="Slope parameter of the logistic function")
    parser.add_argument("-p", '--perform_blasts', dest='PerformBlast', default=True, action="store_false", help="Flag to turn off alignments using blastp. If this is false then alignment directory must be specified") #TODO add ability to specify previously completed alignments.
    parser.add_argument("--blastdir", dest="BlastDir", default=None, help="Directory where the blastp alignments are stored. Must be specified if -p is set.")
    parser.add_argument("-clean", dest="clean", default=False, action='store_true', help="Remove temporary directory. MAKE SURE YOU TRULY WANT TO DO THIS.")
    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument("-i", '--neopred_in', dest="neoPredIn", default=None, type=str,help="Input neoantigen predictions, must be unfiltered or filtered on binding affinity. Example: -I ./AllSamples.neoantigens.txt")
    requiredNamed.add_argument("-f", '--fastas', dest="fastaDir", default="/Users/schencro/Desktop/ChandlerTrevAdInCar/NeoPredPipe/fastaFiles/", type=str,help="Fasta files directory associated with the input.")
    Options = parser.parse_args()  # main user args

    if not Options.neoPredIn:
        parser.error("Some of the required arguments were not provided. Please check required arguments.")
    if len(Options.neorecoOut)!=0:
        if Options.neorecoOut[len(Options.neorecoOut)-1]!='/':
            Options.neorecoOut = Options.neorecoOut+'/'
    if len(Options.fastaDir)!=0:
        if Options.fastaDir[len(Options.fastaDir)-1]!='/':
            Options.fastaDir = Options.fastaDir+'/'
    if Options.PerformBlast==False:
        if Options.BlastDir[len(Options.BlastDir)-1]!='/':
            Options.BlastDir = Options.BlastDir+'/'
        assert Options.BlastDir is not None, "ERROR: A directory must be specified with --blastdir if -p flag is set."
    if Options.clean:
        Clean(Options.neorecoOut)
        sys.exit()

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

def Clean(outputDir):
    os.system('rm -r %s'%(outputDir+"NeoRecoTMP"))

def main():
    print("INFO: Begin.")
    # Pull information about usr system files
    localpath = os.path.abspath(__file__).replace('NeoRecoPo.py', '')  # path to scripts working directory
    Config = configparser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    Options = Parser()
    if Options.Indels:
        tmpFolderName = "NeoRecoTMPIndels"
    else:
        tmpFolderName = "NeoRecoTMP"
    if Options.neorecoOut != "":
        tmpOut = '%s%s/' % (Options.neorecoOut, tmpFolderName)
    else:
        tmpOut = '%s%s%s/' % (Options.neorecoOut, localpath, tmpFolderName)

    netMHCpanPaths = ConfigSectionMap(Config.sections()[1], Config)  # get annovar script paths
    blastp = ConfigSectionMap(Config.sections()[3], Config)['blastp'] # Get blastp from usr_paths file.

    pickleFile = '%s/neorecopo.p'%(tmpOut)
    if os.path.isdir(tmpOut)==False:
        os.system('mkdir %s'%(tmpOut))

    if os.path.isfile(pickleFile) == False:
        preds = StandardPreds(Options) # Create instance of StandardPreds Class
        preds.load() # Load the neoantigen predictions data
        if Options.Indels==False:
            preds.ConstructWTFastas()
            preds.GetWildTypePredictions(netMHCpanPaths) # Extracts the proper WT 'neoantigen'

        preds.BuildFinalTable(Options.Indels)
        with open(pickleFile, 'wb') as outPickle:
            pickle.dump(preds, outPickle)
    else:
        print("INFO: Temporary class file found (neorecopo.p), loading previously processed data.")
        with open(pickleFile,'rb') as inPickle:
            preds = pickle.load(inPickle)
            if preds.wildtypePreds is None and Options.Indels==False:
                preds.GetWildTypePredictions(netMHCpanPaths)
            elif preds.WTandMTtable is None:
                preds.BuildFinalTable(Options.Indels)
            else:
                pass

    if Options.PerformBlast:
        if os.path.isdir(tmpOut):
            if os.path.isdir('%sblastp_results/' % (tmpOut))==False:
                os.system('mkdir %sblastp_results/' % (tmpOut))
        else:
            sys.exit("ERROR: Unable to create blastp_results directory in NeoRecoTMP/")

        preds.PrepBlastPFastaFiles('%sblastp_results/' % (tmpOut))
        preds.PerformBlastPAlignments(blastp, '%sblastp_results/' % (tmpOut))
        preds.PerformCalculations(tmpOut, Options)
        print("INFO: Final table constructed and written.")


    if Options.Dirty:
        os.system('rm -r %s'%(tmpOut))

    print("INFO: Complete.")


if __name__=="__main__":
    main()
