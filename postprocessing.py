#!/usr/bin/env python

import sys
import os
from collections import OrderedDict
import re

def DigestIndSample(toDigest, patName):
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
                        lines.append('\t'.join(line.split()))
                except IndexError as e:
                    pass
    print("INFO: Object size of neoantigens: %s Kb"%(sys.getsizeof(lines)))
    return(lines)

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

    newLines = []
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
        # Step 1: Determine the mutation location in the genotype field
        # Step 2: Output a binary code for present absence in each region for neoantigen heterogeneity
        vcfLine = avReadyLine.split('\t')[8:]
        genoTypeLines = vcfLine[9:] # Extract region specific information


        if Options.colRegions is not None:
            genoTypesPresent = []
            genoTypes = OrderedDict()
            for i in Options.colRegions:
                try:
                    # Only works if the only non-digit characters are base pairs
                    match = re.findall("\:[ACGT]*\:",genoTypeLines[int(i)])[0].replace(":","")
                    if len(match) > 1:
                        present = 1
                    else:
                        present = 0
                    genoTypes.update({'Region_%s'%(i):present})
                    genoTypesPresent.append("+")
                except IndexError as e:
                    genoTypes.update({'Region_%s'%(i):0})
                    genoTypesPresent.append("-")

            regionInfo = '\t'.join([str(genoTypes[i]) for i in genoTypes])
            newLines.append('\t'.join([patName, regionInfo,'line%s' % (epID), chrom, pos, ref, alt, genes, ep]))
        else:
            newLines.append('\t'.join([patName, 'line%s' % (epID), chrom, pos, ref, alt, genes, ep]))

    return(newLines, genoTypesPresent)