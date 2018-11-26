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

def DigestIndSample(toDigest, patName, checkPeptides, pepmatchPaths, indels=False):
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
        pmInputFile = 'tmp/'+patName+'.epitopes.Indels.peptidematch.input'
        pmOutFile = 'tmp/'+patName+'.epitopes.Indels.peptidematch.out'
    else:
        pmInputFile = 'tmp/'+patName+'.epitopes.peptidematch.input'
        pmOutFile = 'tmp/'+patName+'.epitopes.peptidematch.out'

    pmInput = open(pmInputFile,'w')
    for epFile in toDigest:
        print("INFO: Digesting neoantigens for %s" % (patName))
        with open(epFile, 'r') as digest_in:
            for line in digest_in:
                line = line.rstrip('\n')
                try:
                    if line.strip()[0].isdigit():
                        linespl = line.split()
                        lines.append('\t'.join(linespl))
                        if checkPeptides:
                            pmInput.write('>' + linespl[10] + ';' + linespl[2] + '\n' + linespl[2] + '\n')
                except IndexError as e:
                    pass
    pmInput.close()
    if checkPeptides:
        RunPepmatch(pmInputFile, pepmatchPaths['peptidematch_jar'], pepmatchPaths['reference_index'], pmOutFile)
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
    elif 'A' in formatInfo: #alleles found in sample
        genotypeIndex = formatInfo.index('A')
        genotypeFormat = 'allele'
    elif 'AD' in formatInfo: #allele depth info for each allele in sample
        genotypeIndex = formatInfo.index('AD')
        genotypeFormat = 'alldepths'
    else:
        print('INFO: Unknown format in VCF genotype fields, region specific information will probably be incorrect.')
    return(genotypeFormat, genotypeIndex)


def AppendDigestedEps(digestedEps, patName, exonicVars, avReady, Options):
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

    newLines = []
    genoTypesPresent = []
    for ep in digestedEps:
        epID = int(ep.split('\t')[10].split('_')[0].replace('line',''))
        exonicLine = exonicInfo[epID]
        avReadyLine = varInfo[epID]
        chrom = avReadyLine.split('\t')[0]
        pos = avReadyLine.split('\t')[1]
        ref = avReadyLine.split('\t')[3]
        alt = avReadyLine.split('\t')[4]
        genes = ','.join([':'.join(item.split(':')[0:2]) for item in exonicLine.split('\t')[2].split(',') if item != ''])

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
                        match = genoTypeLines[int(i)].split(':')[genotypeIndex]
                        if genotypeFormat == 'allele': #match format: A or AG
                            present = int(len(match) > 1) #present if more than one allele at variant position
                        if genotypeFormat == 'numvarreads': #match format: 0 or 12
                            present = int(int(match) > 0) #present if number of variant reads is > 0
                        if genotypeFormat == 'alldepths': #match format: 10,0 or 10,5
                            present = int(int(match.split(',')[1]) > 0) #present if depth for variant allele > 0
                        genoTypes.update({'Region_%s'%(i):present})
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


def RunPepmatch(pmInput, pepmatchJar, refIndex, pmfileName):
    with open('logForPeptideMatch.tmp', 'a') as logFile:
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
        epkey = line.split('\t')[10]+';'+line.split('\t')[2]
        novel = int(pmDict[epkey]=='No match')
        appendedLines.append(line+'\t'+str(novel))

    return(appendedLines)
