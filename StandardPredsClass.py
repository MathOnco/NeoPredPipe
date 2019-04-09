#!/usr/bin/env python

'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
Contributions from: Eszter Lakatos
'''

import os
import sys
import multiprocessing

from Bio import SeqIO

from NeoAlign import Aligner
from NeoClass import Neoantigen

from vcf_manipulate import ExtractSeq
from predict_binding import predict_neoantigensWT
from postprocessing import DigestIndSampleWT


class StandardPreds:
    '''
    Holds standard predicitons from the other pipeline
    '''
    def __init__(self, Options):
        self.filename = Options.neoPredIn
        self.fastaPath = Options.fastaDir
        self.OutDir = Options.neorecoOut
        self.tmpDir = None
        self.samples = []
        self.hlas = None
        self.fastas = None
        self.filteredPreds = None # Variable containing the predicitons of the mutated sequence
        self.wildtypePreds = None # Variable containing the predictions of the wildtype sequences
        self.tableOutReady = None # Variable containing the WT and MUT sequence table for recognition potential
        self.chop_scores = None # Can be defined by the user by feeding chopscores into the DefinedChopScore function.
        self.WTandMTtable = None # Table containing scores for WT and MT (part of the input of the recognition potential)
        self.ReadyForBlastpFastas = None # List containing fastas for each sample ready to be blasted against IEDB epitopes
        self.blastpResults = None # List containing xml files for the blastp results from each patient epitopes to the IEDB database
        self.filesToPredict = None
        self.epitopeLengths = None
        self.hasExpression = False # Boolean of whether the original predictions had expression information, default is false
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
            if '<=' in line.split('\t'):
                hla = line.split('\t')[-15] #According to netMHCpan output with -BA setting; TODO: make non-hardcoded version
            else:
                hla = line.split('\t')[-13]
            if hla not in self.hlas[sam]:
                self.hlas[sam].append(hla)

        if ':' not in lines[0].split('\t')[-17]: #Checks if it is expression instead of geneName:geneID combination
            self.hasExpression = True


    def __ensureFiltered(self, data):
        '''
        Ensures that neoantigens are <= 500nM binding affinity based on predictions.

        :param data: lines from the input file.
        :return: filtered lines
        '''
        dataOut = []
        has_novelty = data[0].split('\t')[-1] in ['1','0'] #check if pep-match column is added to the end
        for line in data:
            novel = 1
            line = line.split('\t')
            if has_novelty:
                novel = int(line[-1])
                line = line[0:-1]
            if line[len(line)-2]=='<=':
                tmpLine = line[0:len(line)-2]
            else:
                tmpLine = line
            ba = float(tmpLine[len(tmpLine)-2])
            refall = tmpLine[4]
            altall = tmpLine[5]
            if ba <= 500.0 and novel==1:
                dataOut.append('\t'.join(line))

        if os.path.isfile(self.filename.replace('.unfiltered.txt','.filtered.txt'))==False:
            with open(self.filename.replace('.unfiltered.txt','.filtered.txt'), 'w') as outFile:
                outFile.write('\n'.join(dataOut))
        else:
            pass

        return(dataOut)

    def __getEpitopeLengths(self):
        '''
        Gets the epitope lengths for all of the neoantigens within each sample to construct the fasta files.
        :return: Dictionary of samples and the lengths of epitopes that predictions will be made on for the WT.
        '''
        epitopeLengths = {sam: [] for sam in self.samples}
        for k, neo in enumerate(self.filteredPreds):
            neo = neo.split('\t')
            # Unknown number of genotype cols and length may have <= and 'SB'
            if neo[len(neo) - 2] == '<=':
                tmpNeo = neo[0:len(neo) - 2]
            else:
                tmpNeo = neo
            sam = neo[0]
            toMatchFasta = tmpNeo[len(tmpNeo) - 4]
            fasta_head = ';'.join(toMatchFasta.split('_', 1))
            epitope = tmpNeo[len(tmpNeo) - 12]
            epitopeLength = len(tmpNeo[len(tmpNeo) - 12])
            epitopeLengths[sam].append(epitopeLength)

        self.epitopeLengths = epitopeLengths

    def __prepWildtypeFastas(self, tmpDir):
        '''
        Removes and creates new fasta files for the wildtypes. Simple file I/O function.

        :param epitopeLengths: Lengths of epitopes for each sample.
        :param tmpDir: Directory where tempoorary files are stored
        :return: A list of all of the fasta files.
        '''
        # Creates fasta files for wildtype epitopes
        filesToPredict = []
        for sam in self.samples:
            self.epitopeLengths[sam] = list(set(self.epitopeLengths[sam]))
            for epi in self.epitopeLengths[sam]:
                if os.path.isfile('%s%s.wildtype.tmp.%s.fasta' % (tmpDir, sam, epi)):
                    os.system('rm %s%s.wildtype.tmp.%s.fasta' % (tmpDir, sam, epi))
                os.system('touch %s%s.wildtype.tmp.%s.fasta' % (tmpDir, sam, epi))
                filesToPredict.append('%s%s.wildtype.tmp.%s.fasta' % (tmpDir, sam, epi))

        self.filesToPredict = filesToPredict

    def __addToFastaFile(self):
        '''
        Gets the sequence for the wildtype counterpart for the mutant.

        :return:
        '''
        for k, neo in enumerate(self.filteredPreds):
            neo = neo.split('\t')

            # Unknown number of genotype cols and length may have <= and 'SB'
            if neo[len(neo) - 2] == '<=':
                tmpNeo = neo[0:len(neo) - 2]
            else:
                tmpNeo = neo

            sam = neo[0]
            toMatchFasta = tmpNeo[len(tmpNeo) - 4]
            fasta_head = ';'.join(toMatchFasta.split('_', 1))
            epitope = tmpNeo[len(tmpNeo) - 12]
            epitopeLength = len(tmpNeo[len(tmpNeo) - 12])

            seqID, seq = self.__extractSeq(sam, fasta_head, epitopeLength)  # WT seqID and seq

            with open('%s%s.wildtype.tmp.%s.fasta' % (self.tmpDir, sam, epitopeLength), 'a') as tmpFastaOut:
                # tmpFastaOut.write('>' + seqID + sam + ';MT_' + epitope + ";" + str(k) + '\n')
                tmpFastaOut.write('>'+str(k)+';WT'+'\n')
                tmpFastaOut.write(seq + '\n')

        return(None)

    def ConstructWTFastas(self):
        '''
        Orchestrates the creation of the wildtype fasta files.
        :return:
        '''

        self.tmpDir = self.OutDir + 'NeoRecoTMP/'

        if os.path.isdir(self.tmpDir):
            # os.system('rm -r %s' % (tmpDir))
            # os.system('mkdir %s' % (tmpDir))
            pass
        else:
            os.system('mkdir %s' % (self.tmpDir))

        self.__getEpitopeLengths()
        self.__prepWildtypeFastas(self.tmpDir)
        self.__addToFastaFile()



    def GetWildTypePredictions(self, netMHCpan):
        '''
        First it gets the sequence records for both the WT and the MUT AA. It then extracts the proper k-mer from the WT
        for the corresponding mut. Once completed, this is then put into a tmp file for predictions on only the WT
        were the MUT is a predicted neoantigen.
        Holds the following: recordsToGet which is a dictionary with sample : {index of neoantigen : identifier}}

        :netMHCpan: Variable containing the paths loaded from config for netMHCpan
        :return: None, it sets the StandardPreds wildtypePreds variable
        '''

        epcalls = [] # Returns a list of files

        for predictFile in self.filesToPredict:
            # print(predictFile)
            patName = predictFile.split('/')[len(predictFile.split('/'))-1].split('.wildtype.',1)[0]
            hlasNormed = [hla.replace('*','') for hla in self.hlas[patName]]
            epitopeLengths = [predictFile.split('/')[len(predictFile.split('/'))-1].split('.wildtype.')[1].split('.')[1]]
            inFile = {epitopeLengths[0]:predictFile}
            indCalls = predict_neoantigensWT(self.tmpDir, patName, inFile, hlasNormed, epitopeLengths, netMHCpan)
            if indCalls != []:
                indCalls = indCalls[0]
                epcalls.append(indCalls)
            else:
                pass

        filesFromPredictions = {sam:[] for sam in self.samples}
        for rawPreds in epcalls:
            filesFromPredictions[rawPreds.split('/')[len(rawPreds.split('/'))-1].split('.wildtype.',1)[0]].append(rawPreds)

        wildtype_preds = []
        for sample in filesFromPredictions:
            digestedLines = DigestIndSampleWT(filesFromPredictions[sample], sample, False, None)
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
        # test = []
        for pred in self.wildtypePreds:
            item = pred.split('\t')
            wtKey = ','.join([item[0],item[1],item[2],item[11].split('_')[0],str(len(item[3]))])
            outDict.update({wtKey:pred})
            # test.append(wtKey)
        # print(len(test))
        # print(len(list(set(test))))
        return(outDict)

    def BuildFinalTable(self, indels=False):
        '''
        Constructs the final table needed for the recognition potential calculations.

        :return: a list of data needed for the neoantigen recognition potential.
        '''
        # Header for the table is as follows:
        # ID,MUTATION_ID,Sample,WT.PEPTIDE,MT.PEPTIDE,MT.ALLELE (which is the hla for the peptide),WT.SCORE,MT.SCORE,HLA,CHOP_SCORE,Expression (if available)
        if self.chop_scores is None:
            self.SetChopScore(None)
        else:
            pass

        if indels==False:
            wildtypeDict = self.__buildwildtypedict()

        tableLines = []
        count = 1
        k = 0

        for i, MutPred in enumerate(self.filteredPreds):
            MutPred = MutPred.split('\t')

            # Get Mut info
            # Indexes are counted from the end, as number of regions is not fixed
            sample, frame, identifier, mutscore, mutpeptide, hla = MutPred[0], MutPred[-16], MutPred[-6], MutPred[-4], MutPred[-14], MutPred[-15]

            mutKey = ','.join([sample,frame,hla,str(i),str(len(mutpeptide))])

            if indels==False:
                # Get wildtype info
                wtPred = wildtypeDict[mutKey]
                # Check if HLAs match and that the peptide only has one difference in AA
                assert wtPred.split('\t')[2]==hla,"ERROR: HLA types do not match."
                alignDiff = [1 for i in range(0,len(mutpeptide)) if mutpeptide[i] != wtPred.split('\t')[3][i]]
                assert len(alignDiff)==1,"ERROR: Mismatched peptides between Wild Type and Mutant."

                wtscore, wtpeptide = wtPred.split('\t')[13], wtPred.split('\t')[3]
            else:
                wtscore = "1000.0"
                wtpeptide = "-"

            # Get the HLAs for this sample
            lineHLAs = '"' + ','.join([sampleHLA.replace("HLA-","").replace("*","").replace(":","") for sampleHLA in self.hlas[sample]]) + '"'

            # Construct table line
            lineOut = '\t'.join([str(count), identifier, sample, wtpeptide, mutpeptide, hla.replace("HLA-","").replace("*","").replace(":",""), wtscore, mutscore, lineHLAs, str(self.chop_scores[count-1])])
            if self.hasExpression:
                expression = MutPred[-17]
                lineOut = '\t'.join([lineOut, expression])
            count += 1
            tableLines.append(lineOut)

        if indels:
            print("INFO: Wild-type affinity values for indels are fixed to %s."%wtscore)
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

    def PerformBlastPAlignments(self, blastp, outputDir):
        '''
        Performs alignments from IEDB sequences and predicted neoantigens. It constructs xml files that are stored in
        the temporary directory.

        :param blastp: usr_paths.ini executable for NCBI's blastp
        :param outputDir: temporary directory where blastp xml result files are stored.
        :return: None. It sets self.blastpResults to a list of xml files for each patient.
        '''

        blastpOut = []
        for fasta in self.ReadyForBlastpFastas:
            assert os.path.isfile(fasta), "ERROR: Unable to locate fasta file %s"%(fasta)
            out = "%s.blastpResults.xml"%(fasta.replace('.readyForBlastp.fasta',''))
            iedb = os.path.realpath(__file__).split('StandardPredsClass')[0] + "ncbi_epitope_db/IEDB_positive_T-cell_assays.fasta"

            if os.path.isfile(out)==False:
                print("INFO: Running blastp on %s"%(fasta))
                cmd = ' '.join([blastp, '-evalue 100000000 -gapextend 1 -gapopen 11 -outfmt 5 -out', out, '-query', fasta, '-subject', iedb])
                os.system(cmd)
                blastpOut.append(out)
            else:
                print("INFO: Found blastp results %s. Nothing to be done."%(fasta))

        self.blastpResults = blastpOut

    def PrepBlastPFastaFiles(self, outputDir):
        '''
        Constructs a fasta file for each patients neoantigens MUT and corresponding WT epitopes to feed into blast.

        :param outputDir: Directory that houses the NeoReco information.
        :return: None. It sets self.ReadyForBlastpFastas equal to the  Fasta file for each of the patients.
        '''
        sample_epitopes = dict.fromkeys([item.split('\t')[2] for item in self.WTandMTtable])
        for sample in sample_epitopes:
            sample_epitopes[sample]=[]

        for entry in self.WTandMTtable:
            sample_epitopes[entry.split('\t')[2]].append(entry)

        fastaFiles = []
        for sample in sample_epitopes:
            outputFile = '%s%s.readyForBlastp.fasta'%(outputDir,sample)
            fastaFiles.append(outputFile)
            if os.path.isfile(outputFile) == False:
                with open(outputFile, 'w') as usrOut:
                    for entry in sample_epitopes[sample]:
                        entry = entry.split('\t')
                        usrOut.write('|'.join(['>%s'%(sample),entry[1],'MT',entry[0]]) + '\n')
                        usrOut.write(entry[4] + '\n')
            else:
                pass

        self.ReadyForBlastpFastas = fastaFiles

    def _buildNeoFile(self, tmpDir):
        '''
        Writes the neoantigen table for the calculations of the neoantigen recognition potential. Writes the information to the tmp directory.

        :param tmpDir: Directory housing the temporary files.
        :return: None
        '''
        filename = "%sNeoantigens.WTandMTtable.txt"%(tmpDir)

        with open(filename, 'w') as outFile:
            if len(self.WTandMTtable[0].split('\t'))==11:
                outFile.write("\t".join(["ID","MUTATION_ID","Sample","WT.PEPTIDE","MT.PEPTIDE","MT.ALLELE","WT.SCORE","MT.SCORE","HLA","CHOP_SCORE", "EXPRESSION"]) + "\n")
            else:
                outFile.write("\t".join(["ID","MUTATION_ID","Sample","WT.PEPTIDE","MT.PEPTIDE","MT.ALLELE","WT.SCORE","MT.SCORE","HLA","CHOP_SCORE"]) + "\n")
            for line in self.WTandMTtable:
                outFile.write(line + "\n")

        return(filename)

    def _compileNeoantigens(self, neofile, indels):
        '''
        Reads and builds NeoClass class

        :param neofile: File constructed with _buildNeoFile
        :return: neoantigens, samples
        '''
        neoantigens = {}
        with open(neofile, 'r') as f:
            header = f.readline()
            htab = header.strip().split("\t")
            hdict = {}
            for i in range(0, len(htab)):
                hdict[htab[i]] = i

            line = f.readline()
            while line:
                line = line.strip()
                nparams = line.split("\t")
                if nparams[7] == "NA":
                    line = f.readline()
                    continue
                neoantigen = Neoantigen(nparams, indels)

                neoantigens[neoantigen.id] = neoantigen
                neoantigen.setA()
                line = f.readline()

        samples = set(map(lambda neo: neo.getSampleName(), neoantigens.values()))
        return [neoantigens, samples]

    def PerformCalculations(self, tmpDir, Options):
        '''
        Main orchestration for performing the neoantigen recognition potential.

        :param tmpDir: Directory of the temporary files
        :param Options: Config class holding the outputdir and parameters for a and k
        :return: None. Writes the final output.
        '''

        a = Options.a # Midpoint parameter of the logistic function, alignment score threshold
        k = Options.k # Slope parameter of the logistic function
        neofile = self._buildNeoFile(tmpDir)

        outFile = Options.neorecoOut # Will output to working directory if no output directory specified

        [neoantigens, samples] = self._compileNeoantigens(neofile, Options.Indels)

        aligner = Aligner()

        for sample in samples:
            xmlpath = "%sblastp_results/%s.blastpResults.xml"%(tmpDir,sample)
            assert os.path.isfile(xmlpath), "ERROR: BlastP xml file %s not found"%(xmlpath)

            aligner.readAllBlastAlignments(xmlpath)

        aligner.computeR(a, k)

        if Options.Indels:
            outputFile = "PredictedRecognitionPotentials.Indels"
        else:
            outputFile = "PredictedRecognitionPotentials"

        with open(Options.neorecoOut + outputFile + ".txt", "w") as outFile:
            if neoantigens[neoantigens.keys()[0]].expr is None:
                header = ["NeoantigenID", "Mutation", "Sample", "MutatedPeptide", "ResidueChangeClass", "MutantPeptide",
                "WildtypePeptide", "HLA", "A", "R", "Excluded", "NeoantigenRecognitionPotential"]
            else:
                header = ["NeoantigenID", "Mutation", "Sample", "MutatedPeptide", "ResidueChangeClass", "MutantPeptide",
                "WildtypePeptide", "HLA","Expression", "A", "R", "Excluded", "NeoantigenRecognitionPotential"]
            header = "\t".join(header)
            outFile.write(header+'\n')
            for neo in neoantigens:
                neoantigen = neoantigens[neo]
                w = neoantigen.getWeight()  # excludes neoantigens that mutated from a nonhydrophobic residue on position 2 or 9
                A = neoantigens[neo].getA()  # MHC amplitude A
                mtpeptide = neoantigens[neo].mtPeptide  # mutant peptide
                wtpeptide = neoantigens[neo].wtPeptide
                R = aligner.getR(neo)

                # Residue change:
                # HH: from hydrophobic to hydrophobic,
                # NN: from non-hydrophobic to non-hydrophobic
                # HN: from hydrophobic to non-hydrophobic,
                # NH: from non-hydrophobic to hydrophobic
                # other (WW, WH, HW, NW, WN) which include aminoacids without a clear classification
                residueChange = neoantigen.residueChange

                fitnessCost = A * R * w

                if neoantigen.expr is not None:
                    l = [neo, neoantigen.mid, neoantigen.sample, neoantigen.position, residueChange, mtpeptide, wtpeptide, neoantigen.allele,
                    neoantigen.expr, A,R, 1 - w, fitnessCost]
                else:
                    l = [neo, neoantigen.mid, neoantigen.sample, neoantigen.position, residueChange, mtpeptide, wtpeptide, neoantigen.allele, A,
                    R, 1 - w, fitnessCost]  
                l = "\t".join(map(lambda s: str(s), l))
                outFile.write(l+'\n')
