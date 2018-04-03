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
from vcf_manipulate import ExtractSeq, predict_neoantigensWT
from postprocessing import DigestIndSample

def Parser():
    # get user variables
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dirty", dest="Dirty", default=True, action='store_false', help="Flag to keep intermediate files. Default: False. Note: Use for debugging.")
    parser.add_argument("-o", '--neoreco_out', dest="neorecoOut", default="", type=str, help="Output Directory.")
    parser.add_argument("-a", '--midpoint', dest='a', default=1., type=float, help="Midpoint parameter of the logistic function, alignment score threshold.")
    parser.add_argument("-k", '--slope', dest='k', default=1., type=float, help="Slope parameter of the logistic function")
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

class StandardPreds:
    '''
    Holds standard predicitons from the other pipeline
    '''
    def __init__(self, Options):
        self.filename = Options.neoPredIn
        self.fastaPath = Options.fastaDir
        self.OutDir = Options.neorecoOut
        self.samples = []
        self.hlas = None
        self.fastas = None
        self.filteredPreds = None # Variable containing the predicitons of the mutated sequence
        self.wildtypePreds = None # Variable containing the predictions of the wildtype sequences
        self.tableOutReady = None # Variable containing the WT and MUT sequence table for recognition potential
        self.chop_scores = None # Can be defined by the user by feeding chopscores into the DefinedChopScore function.
        self.WTandMTtable = None # Table containing scores for WT and MT (part of the input of the recognition potential)

    def load(self):
        '''
        Loads the data class of neoantigen predictions to get the right information.

        :return: None
        '''
        with open(self.filename, 'r') as inFile:
            lines = [line.replace('\n','') for line in inFile.readlines()]

        lines = self.__ensureFiltered(lines)

        self.filteredPreds = lines
        self.samples = list(set([line.split('\t')[0] for line in lines]))
        self.hlas = {sam:[] for sam in self.samples}
        self.fastas = {sam:'%s%s.reformat.fasta'%(self.fastaPath,sam) for sam in self.samples}

        for line in lines:
            sam = line.split('\t')[0]
            hla = line.split('\t')[11]
            if hla not in self.hlas[sam]:
                self.hlas[sam].append(hla)

    def __ensureFiltered(self, data):
        '''
        Ensures that neoantigens are <= 500nM binding affinity based on predictions.

        :param data: lines from the input file.
        :return: filtered lines
        '''
        dataOut = []
        for line in data:
            line = line.split('\t')
            if line[len(line)-2]=='<=':
                tmpLine = line[0:len(line)-2]
            else:
                tmpLine = line
            ba = float(tmpLine[len(tmpLine)-2])
            if ba <= 500.0:
                dataOut.append('\t'.join(line))
        return(dataOut)

    def GetWildTypePredictions(self, netMHCpan):
        '''
        First it gets the sequence records for both the WT and the MUT AA. It then extracts the proper k-mer from the WT
        for the corresponding mut. Once completed, this is then put into a tmp file for predictions on only the WT
        were the MUT is a predicted neoantigen.
        Holds the following: recordsToGet which is a dictionary with sample : {index of neoantigen : identifier}}

        :netMHCpan: Variable containing the paths loaded from config for netMHCpan
        :return: None, it sets the StandardPreds wildtypePreds variable
        '''
        tmpDir = self.OutDir + 'NeoRecoTMP/'
        if os.path.isdir(tmpDir):
            os.system('rm -r %s' % (tmpDir))
            os.system('mkdir %s' % (tmpDir))
        else:
            os.system('mkdir %s' % (tmpDir))

        epitopeLengths = {sam:[] for sam in self.samples}
        for neo in self.filteredPreds:
            neo = neo.split('\t')
            # Unknown number of genotype cols and length may have <= and 'SB'
            if neo[len(neo)-2]=='<=':
                tmpNeo = neo[0:len(neo)-2]
            else:
                tmpNeo = neo
            sam = neo[0]
            toMatchFasta = tmpNeo[len(tmpNeo)-4]
            fasta_head = ';'.join(toMatchFasta.split('_',1))
            epitope = tmpNeo[len(tmpNeo)-5]
            epitopeLength = len(tmpNeo[len(tmpNeo)-5])
            epitopeLengths[sam].append(epitopeLength)

        # Creates fasta files for wildtype epitopes
        filesToPredict = []
        for sam in self.samples:
            epitopeLengths[sam] = list(set(epitopeLengths[sam]))
            for epi in epitopeLengths[sam]:
                if os.path.isfile('%s%s.wildtype.tmp.%s.fasta' % (tmpDir, sam, epi)):
                    os.system('rm %s%s.wildtype.tmp.%s.fasta' % (tmpDir, sam, epi))
                os.system('touch %s%s.wildtype.tmp.%s.fasta' % (tmpDir, sam, epi))
                filesToPredict.append('%s%s.wildtype.tmp.%s.fasta' % (tmpDir, sam, epi))

        seen = []
        for neo in self.filteredPreds:
            neo = neo.split('\t')
            # Unknown number of genotype cols and length may have <= and 'SB'
            if neo[len(neo)-2]=='<=':
                tmpNeo = neo[0:len(neo)-2]
            else:
                tmpNeo = neo
            sam = neo[0]
            toMatchFasta = tmpNeo[len(tmpNeo)-4]
            fasta_head = ';'.join(toMatchFasta.split('_',1))
            epitope = tmpNeo[len(tmpNeo)-5]
            epitopeLength = len(tmpNeo[len(tmpNeo)-5])
            seqID, seq = self.__extractSeq(sam, fasta_head, epitopeLength) # WT seqID and seq

            if seqID+seq+str(epitopeLength) in seen:
                pass
            else:
                seen.append(seqID+seq+str(epitopeLength))
                with open('%s%s.wildtype.tmp.%s.fasta'%(tmpDir,sam,epitopeLength), 'a') as tmpFastaOut:
                    tmpFastaOut.write('>' + seqID + '\n')
                    tmpFastaOut.write(seq + '\n')

        epcalls = [] # Returns a list of files
        for predictFile in filesToPredict:
            patName = predictFile.split('/')[len(predictFile.split('/'))-1].split('.',1)[0]
            hlasNormed = [hla.replace('*','') for hla in self.hlas[patName]]
            epitopeLengths = [predictFile.split('/')[len(predictFile.split('/'))-1].split('.')[3]]
            inFile = {epitopeLengths[0]:predictFile}
            epcalls.append(predict_neoantigensWT(tmpDir, patName, inFile, hlasNormed, epitopeLengths, netMHCpan)[0])

        filesFromPredictions = {sam:[] for sam in self.samples}
        for rawPreds in epcalls:
            filesFromPredictions[rawPreds.split('/')[len(rawPreds.split('/'))-1].split('.',1)[0]].append(rawPreds)

        wildtype_preds = []
        for sample in filesFromPredictions:
            digestedLines = DigestIndSample(filesFromPredictions[sample], sample, False, None)
            appendedLines = []
            for line in digestedLines:
                appendedLines.append('\t'.join([sample,line]))
            wildtype_preds = wildtype_preds + appendedLines

        self.wildtypePreds = wildtype_preds

        return(None)

    def __extractSeq(self, sample, identifier, epitopeLength):
        '''
        Extracts the sequence from the *.tmp.epi.fasta file and reverts the sequence back.

        :return: the wildtype sequence and header
        '''
        WT = []
        Mut = []
        count=0
        for seq_record in SeqIO.parse(self.fastas[sample], 'fasta'):
            seqIdentifier = ';'.join(seq_record.id.split(';',3)[0:2])[0:16]
            if identifier in seqIdentifier:
                count += 1

                if 'WILDTYPE' in seq_record.id.split(';')[2]:
                    WT.append(seq_record.id)
                    WT.append(seq_record)
                else:
                    try:
                        pos = int(seq_record.id.replace(";;", ";").split(";")[5]) - 1
                    except ValueError:
                        pos = int(seq_record.id.replace(";;", ";").split(";")[6]) - 1
                    Mut.append(seq_record.id)
                    Mut.append(seq_record)

                if count==2:
                    break

        WTepiSeq = ExtractSeq(WT[1], pos, epitopeLength)

        return(WT[0], WTepiSeq)

    def __buildwildtypedict(self):
        '''
        Constructs a dictionary for wildtype predictions to be extracted from to check for a match

        :return: a dictionary of wildtype predictions to extract from
        '''
        outDict = {}
        for pred in self.wildtypePreds:
            item = pred.split('\t')
            wtKey = ','.join([item[0],item[1],item[2],item[11],str(len(item[3]))])
            outDict.update({wtKey:pred})
        return(outDict)

    def BuildFinalTable(self):
        '''
        Constructs the final table needed for the recognition potential calculations.

        :return: a list of data needed for the neoantigen recognition potential.
        '''
        # Header for the table is as follows:
        # ID,MUTATION_ID,Sample,WT.PEPTIDE,MT.PEPTIDE,MT.ALLELE (which is the hla for the peptide),WT.SCORE,MT.SCORE,HLA,CHOP_SCORE
        if self.chop_scores is None:
            self.SetChopScore(None)
        else:
            pass

        wildtypeDict = self.__buildwildtypedict()

        tableLines = []
        count = 1
        for MutPred in self.filteredPreds:
            MutPred = MutPred.split('\t')

            # Get Mut info
            sample, frame, identifier, mutscore, mutpeptide, hla = MutPred[0], MutPred[10], MutPred[20], MutPred[22], MutPred[12], MutPred[11]

            mutKey = ','.join([sample,frame,hla,identifier,str(len(mutpeptide))])

            # Get wildtype info
            wtPred = wildtypeDict[mutKey]

            # Check if HLAs match and that the peptide only has one difference in AA
            assert wtPred.split('\t')[2]==hla,"ERROR: HLA types do not match."
            alignDiff = [1 for i in range(0,len(mutpeptide)) if mutpeptide[i] != wtPred.split('\t')[3][i]]
            assert len(alignDiff)==1,"ERROR: Mismatched peptides between Wild Type and Mutant."

            wtscore, wtpeptide = wtPred.split('\t')[13], wtPred.split('\t')[3]

            # Get the HLAs for this sample
            lineHLAs = '"' + ','.join([sampleHLA.replace("HLA-","").replace("*","").replace(":","") for sampleHLA in self.hlas[sample]]) + '"'

            # Construct table line
            lineOut = '\t'.join([str(count), identifier, sample, wtpeptide, mutpeptide, hla.replace("HLA-","").replace("*","").replace(":",""), wtscore, mutscore, lineHLAs, str(self.chop_scores[count-1])])
            count += 1
        tableLines.append(lineOut)

        self.WTandMTtable = tableLines


    def SetChopScore(self, chopscores):
        '''
        Adds chopscores for the final input if so desired. Data structure must be a list that matches the indices of StandardPreds.filteredPreds

        :param chopscores: List of integer or floats that match the indices of StandardPreds.filteredPreds
        :return: None, it sets StandardPreds.chop_scores
        '''
        if chopscores is None:
            self.chop_scores = [1 for i in range(0,len(self.filteredPreds))]
        else:
            self.chop_scores=chopscores

def main():
    print("INFO: Begin.")
    # Pull information about usr system files
    localpath = os.path.abspath(__file__).replace('NeoRecoPo.py', '')  # path to scripts working directory
    Config = configparser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    Options = Parser()
    netMHCpanPaths = ConfigSectionMap(Config.sections()[1], Config)  # get annovar script paths

    tmpOut = '%s/%s/NeoRecoTMP/'%(Options.neorecoOut, localpath)
    pickleFile = '%s/neorecopo.p'%(tmpOut)
    if os.path.isdir(tmpOut)==False:
        os.system('mkdir %s'%(tmpOut))

    if os.path.isfile(pickleFile) == False:
        preds = StandardPreds(Options) # Create instance of StandardPreds Class
        preds.load() # Load the neoantigen predictions data
        preds.GetWildTypePredictions(netMHCpanPaths) # Extracts the proper WT 'neoantigen'
        preds.BuildFinalTable()
        with open(pickleFile, 'wb') as outPickle:
            pickle.dump(preds, outPickle)
    else:
        with open(pickleFile,'rb') as inPickle:
            preds = pickle.load(inPickle)
            if preds.wildtypePreds is None:
                preds.GetWildTypePredictions(netMHCpanPaths)
            elif preds.WTandMTtable is None:
                preds.BuildFinalTable()
            else:
                pass


    if Options.Dirty:
        os.system('rm -r %s'%(tmpOut))

    print("INFO: Complete.")


if __name__=="__main__":
    main()