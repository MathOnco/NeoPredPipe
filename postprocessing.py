#!/usr/bin/env python

'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
Contributions from: Eszter Lakatos
'''

import sys
import os
from collections import OrderedDict
import subprocess
import re
from process_expression import BuildExpTable, BuildGeneIDTable

def DigestIndSample(toDigest, patName, Options, pepmatchPaths, indels=False):
    '''
    Filters the resulting file and strips all information within it down to individual calls.

    :param toDigest: A list of files to be digested for an individual patient.
    :param patName: Patient/sample identifier
    :return: All Neoantigen Prediction lines free of other information in prediction files.
    '''
    # temp_files = None
    # output_file = "%s%s.digested.txt" % (FilePath, toDigest[0].split('/')[len(toDigest[0].split('/')) - 1].split('.epitopes.')[0])

    lines = []
    if indels:
        pmInputFile = Options.OutputDir+'tmp/'+patName+'.epitopes.Indels.peptidematch.input'
        pmOutFile = Options.OutputDir+'tmp/'+patName+'.epitopes.Indels.peptidematch.out'
    else:
        pmInputFile = Options.OutputDir+'tmp/'+patName+'.epitopes.peptidematch.input'
        pmOutFile = Options.OutputDir+'tmp/'+patName+'.epitopes.peptidematch.out'

    pmInputDict = {}
    pmInputLines = 0
    for epFile in toDigest:
        print("INFO: Digesting neoantigens for %s" % (patName))
        with open(epFile, 'r') as digest_in:
            for line in digest_in:
                line = line.rstrip('\n')
                try:
                    if line.strip()[0].isdigit():
                        linespl = line.split()
                        lines.append('\t'.join(linespl))
                        if Options.checkPeptides:
                            pep = linespl[2]
                            pmInputDict['>'+pep] = pep
                            pmInputLines+=1
                except IndexError as e:
                    pass
                    
    if Options.checkPeptides and pmInputLines>0:
        RunPepmatch(Options.OutputDir,pmInputDict,pmInputFile, pepmatchPaths['peptidematch_jar'], pepmatchPaths['reference_index'], pmOutFile)
        lines = ProcessPepmatch(pmOutFile, lines)
    print("INFO: Object size of neoantigens: %s Kb"%(sys.getsizeof(lines)))
    return(lines)

def DigestIndSampleWT(toDigest, patName, checkPeptides, pepmatchPaths):
    '''
    Filters the resulting file and strips all information within it down to individual calls.

    :param toDigest: A list of files to be digested for an individual patient.
    :param patName: Patient/sample identifier
    :return: All Neoantigen Prediction lines free of other information in prediction files.
    '''
    # temp_files = None
    # output_file = "%s%s.digested.txt" % (FilePath, toDigest[0].split('/')[len(toDigest[0].split('/')) - 1].split('.epitopes.')[0])

    lines = []
    for epFile in toDigest:
        print("INFO: Digesting neoantigens for %s" % (patName))
        with open(epFile, 'r') as digest_in:
            for line in digest_in:
                line = line.rstrip('\n')
                try:
                    if line.strip()[0].isdigit():
                        linespl = line.split()
                        lines.append('\t'.join(linespl))
                except IndexError as e:
                    pass
    print("INFO: Object size of neoantigens: %s Kb"%(sys.getsizeof(lines)))
    return(lines)

def DefineGenotypeFormat(testLine):
    '''
    Determines which element of genotype fields contains relevant information and in what format
    Current options are NV (number of reads with variant allele) and A (alleles found in sample)

    :param testLine: A single line from exonic_variant_function file
    :return: Genotype format (allele or numvarreads or alldepths) and the corresponding information's index in genotype info
    '''
    genotypeFormat = 'unknown'
    genotypeIndex = -1
    formatInfo = testLine.split('\t')[19].split(':')
    if 'NV' in formatInfo: #number of variant reads in sample
        genotypeIndex = formatInfo.index('NV')
        genotypeFormat = 'numvarreads'
        print('Processing genotype information according to NV (number of variant reads) field.')
    elif 'A' in formatInfo: #alleles found in sample
        genotypeIndex = formatInfo.index('A')
        genotypeFormat = 'allele'
        print('Processing genotype information according to A (list of alleles) field.')
    elif 'FREQ' in formatInfo: #varscan2 allele frequency
        genotypeIndex = formatInfo.index('FREQ')
        genotypeFormat = 'varscanfreq'
        print('Processing genotype information according to FREQ (allele frequency) field.')
    elif 'AD' in formatInfo: #allele depth info for each allele in sample
        genotypeIndex = formatInfo.index('AD')
        genotypeFormat = 'alldepths'
        print('Processing genotype information according to AD (allelic depth) field.')
    elif 'GT' in formatInfo: #just use genotype-info, should be very universal
        genotypeIndex = formatInfo.index('GT')
        genotypeFormat = 'genotype'
        print('Processing genotype information according to GT (genotype) field.')
    elif 'AU' in formatInfo and 'CU' in formatInfo and 'GU' in formatInfo and 'TU' in formatInfo:
        genotypeIndex = {a:formatInfo.index(a) for a in ['AU','CU','GU','TU']} #tier count information for each possible base
        genotypeFormat = 'strelka'
        print('INFO: Strelka-specific format found, processing genotype information accordingly.')
    else:
        print('INFO: Unknown format in VCF genotype fields, region specific information will be missing. See readme for supported formats.')
    return(genotypeFormat, genotypeIndex)

def getAltIndex(allele, alleleList):
    ind = -1
    if allele in alleleList:
        ind = alleleList.index(allele)
    return(ind)


def AppendDigestedEps(FilePath,digestedEps, patName, exonicVars, avReady, Options):
    '''
    Appends information to digested Eps for analysis.

    :param digestedEps: Lines from DigestIndSample from netMHCpan
    :param patName: Patient/sample identifier
    :param exonicVars: exonic_variant_function files from ANNOVAR
    :param avReady: ANNOVAR ready file created from VCF
    :param Options: ArgParse options for processing multiregion
    :return: List of new lines with appended information
    '''

    # Pull exonic_variant_function files into dictionary with line # as key
    with open(exonicVars, 'r') as exInfo:
        lines = exInfo.readlines()
        exonicInfo = {int(line.split('\t')[0].replace('line','')):line.rstrip("\n") for line in lines}
    with open(avReady, 'r') as avIn:
        lines = avIn.readlines()
        varInfo = {i+1:line.rstrip("\n") for i,line in enumerate(lines)}

    # Test one line to determine FORMAT of genotype fields in vcf
    testLine = exonicInfo[min(exonicInfo.keys())]
    if Options.colRegions is not None:
        genotypeFormat, genotypeIndex = DefineGenotypeFormat(testLine)

    # Load in Expression data if available
    expTable = None
    if Options.expression is not None:
        # Options.allExpFiles is either a single string or a list of filename strings
        if isinstance(Options.allExpFiles, list):
            sampleExpFile = filter(lambda x: '/'+patName+'.tsv' in x, Options.allExpFiles)
        else:
            sampleExpFile = [Options.allExpFiles]
        if len(sampleExpFile)>0:
            idType, expTable = BuildExpTable(sampleExpFile[0], Options.expMultiregion) #take zeroth element to unlist, there should not be ambiguity
            idTable = BuildGeneIDTable(FilePath,idType)
        else:
            print('INFO: No expression file found for sample %s!'%(patName))

    newLines = []
    genoTypesPresent = []
    for ep in digestedEps:
        if Options.typeII:
            epID = int(ep.split('\t')[3].split(';')[0].replace('line',''))
        else:
            epID = int(ep.split('\t')[10].split('_')[0].replace('line',''))
        exonicLine = exonicInfo[epID]
        avReadyLine = varInfo[epID]
        chrom = avReadyLine.split('\t')[0]
        pos = avReadyLine.split('\t')[1]
        ref = avReadyLine.split('\t')[3]
        alt = avReadyLine.split('\t')[4]
        altOriginal = avReadyLine.split('\t')[12] # retrieve original alternative allele to handle multi-allele variants
        altOriginal = altOriginal.split(',')
        if len(altOriginal) > 1:
            altInd = getAltIndex(alt,altOriginal) # get which of multiple alleles is this, in case needed
        else:
            altInd = 0
        geneList = [':'.join(item.split(':')[0:2]) for item in exonicLine.split('\t')[2].split(',') if item != '']
        nmList = [item.split(':')[1] for item in geneList]
        genes = ','.join(geneList)
        
        if Options.expression is not None:
            geneExp = 'NA'
            if expTable is not None:
                for nmId in nmList:
                    #convert gene ID in NeoPredPipe to geneID that is found in the exp file and query expression dictionary
                    if nmId in idTable.keys():
                        tableID = idTable[nmId]
                        if tableID in expTable.keys():
                            geneExp = expTable[tableID]
                            break
            genes = '\t'.join([genes, geneExp])

        # Getting information from the genotype fields
        # Step 1: Determine the mutation location in the genotype field (based on FORMAT info/ genotype index)
        # Step 2: Output a binary code for present absence in each region for neoantigen heterogeneity
        vcfLine = avReadyLine.split('\t')[8:]
        genoTypeLines = vcfLine[9:] # Extract region specific information

        if Options.colRegions is not None:
            genoTypesPresent = []
            genoTypes = OrderedDict()
            for i in Options.colRegions:
                if genotypeFormat != 'unknown' and len(genoTypeLines)>int(i):
                    #first handle multiregion vcf where absence is indicated by .
                    if genoTypeLines[int(i)]=='.':
                        genoTypes.update({'Region_%s' % (i): 0})
                    else:
                    #otherwise take the corresponding piece of the genotype line
                        # first handle Strelka, as in this case genotypeIndex is a dictionary, not an integer
                        if genotypeFormat=='strelka': #match format is tier1,tier2 for all 4 alleles
                            match = genoTypeLines[int(i)].split(':')[genotypeIndex[alt+'U']] #use alt to identify which allele to look for
                            #match format: tier1,tier2, so 0,0 or 5,6
                            present = int(sum([int(x) for x in match.split(',')]) > 0) #present if any is above zero, meaning sum is above 0
                            #present = int(match.split(',')[0]) > 0 #present if tier1 is above zero
                        else:
                            match = genoTypeLines[int(i)].split(':')[genotypeIndex]
                        try:
                            if genotypeFormat == 'allele': #match format: A or AG
                                present = int(len(match) > 1) #present if more than one allele at variant position
                            if genotypeFormat == 'numvarreads': #match format: 0 or 12
                                match = [int(x) for x in match.split(',')]
                                if altInd > -1 and len(match)>1:
                                    present = int(match[altInd] > 0) #present if number of variant reads is > 0
                                else:
                                    present = int(max(match) > 0)
                            if genotypeFormat == 'alldepths': #match format: 10,0 or 10,5
                                matchAltStr = [x for x in match.split(',')[1:]] #discard reference allele depth
                                if len(matchAltStr)==0: # handle case when there is only a single '.' listed
                                    present = 0
                                else:
                                    matchAlt = []
                                    for x in matchAltStr:
                                        if x == '.': # handle case when '.' is found instead of 0
                                            matchAlt.append(0)
                                        else:
                                            matchAlt.append(int(x))
                                    if altInd > -1: #might have to handle multi-allele case by finding which of the alleles
                                        present = int(matchAlt[altInd] > 0) #present if depth for variant allele > 0
                                    else:
                                        present = int(max(matchAlt) > 0)
                            if genotypeFormat=='varscanfreq': #match format 18.1% or 0.0%
                                present = int(float(match.replace('%',''))>0) #present if float before % is > 0
                            if genotypeFormat=='genotype': #match format 0/0, 0/1 or 1/0
                                present = int('1' in match) #present if contains 1
                            genoTypes.update({'Region_%s'%(i):present})
                        except IndexError:
                            sys.exit("ERROR: Index problem while Processing: %s" % (avReadyLine))
                    genoTypesPresent.append("+")
                else:
                    genoTypes.update({'Region_%s'%(i):-1})
                    genoTypesPresent.append("-")

            regionInfo = '\t'.join([str(genoTypes[i]) for i in genoTypes])
            newLines.append('\t'.join([patName, regionInfo,'line%s' % (epID), chrom, pos, ref, alt, genes, ep]))
        else:
            genoTypesPresent = []
            newLines.append('\t'.join([patName, 'line%s' % (epID), chrom, pos, ref, alt, genes, ep]))

    return(newLines, genoTypesPresent)


def RunPepmatch(FilePath,pmInputDict, pmInput, pepmatchJar, refIndex, pmfileName):
    with open(pmInput, 'w') as inputFile:
        for k in pmInputDict.keys():
            inputFile.write(k+'\n'+pmInputDict[k]+'\n')

    with open(FilePath+'logForPeptideMatch.tmp', 'a') as logFile:
        cmd = ['java', '-jar', pepmatchJar, '-a', 'query', '-i', refIndex,'-Q', pmInput, '-o', pmfileName]
        runcmd = subprocess.Popen(cmd, stdout=logFile)
        runcmd.wait()

def ProcessPepmatch(pmfileName, epLines):
    with open(pmfileName, 'r') as pmFile:
        pmFile.readline()
        pmFile.readline() #read first two header lines
        pmDict = {line.split('\t')[0] : line.split('\t')[1].rstrip('\n') for line in pmFile.readlines() }
    appendedLines = []
    for line in epLines:
        epkey = line.split('\t')[2]
        novel = int(pmDict[epkey]=='No match')
        appendedLines.append(line+'\t'+str(novel))

    return(appendedLines)
