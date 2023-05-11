#!/usr/bin/env python

'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
Contributions from: Eszter Lakatos

Adapted from Marta Luksza (see RecognitionPotential.md)
'''

from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from math import log, exp
import sys

class Aligner(object):
    '''
    Class to compute alignment score of neoantigens with IEDB epitopes and compute TCR-recognition
    probabilities.
    '''
    INF = float("inf")

    @staticmethod
    def align(seq1, seq2):
        '''
        Smith-Waterman alignment with default parameters.
        '''
        matrix = matlist.blosum62
        gap_open = -11
        gap_extend = -1
        aln = pairwise2.align.localds(seq1.upper(), seq2.upper(), matrix, gap_open, gap_extend)
        return aln

    @staticmethod
    def logSum(v):
        '''
        compute the logarithm of a sum of exponentials

        :param v: binding energies
        :return: log(sum(map(lambda x: exp(x - ma), v))) + ma
        '''

        if len(v) == 0:
            return -Aligner.INF
        ma = max(v)
        if ma == -Aligner.INF:
            return -Aligner.INF
        return log(sum(map(lambda x: exp(x - ma), v))) + ma

    def __init__(self):
        # dictionary of computed Ri-values mapped to neoantigen identifiers
        self.Ri = {}
        # dictionary of IEDB epitope alignments mapped to neoantigen identifiers
        self.alignments = {}
        # dictionary of the highest scoring alignments mapped to neoantigen identifiers
        self.maximum_alignment = {}

    def readAllBlastAlignments(self, xmlpath):
        '''
        Read precomputed blastp alignments from xml files,
        compute alignment scores,
        find the highest scoring alignment for each neoantigen.

        :param xmlpath: File containing pre-blasted neoantigens against the IEDB eptiopes
        :return: None. Sets self.alignments
        '''
        f = open(xmlpath)
        blast_records = NCBIXML.parse(f)

        maxscore = {}
        for brecord in blast_records:
            tab = str(brecord.query).split("|")
            ptype = tab[2]
            nid = int(tab[3])
            if ptype == "MT":
                if not nid in maxscore:
                    maxscore[nid] = 0
                for alignment in brecord.alignments:
                    if not nid in self.alignments:
                        self.alignments[nid] = {}
                        self.maximum_alignment[nid] = None
                        self.maximum_alignment[nid] = 0
                        maxscore[nid] = 0
                    species = " ".join((str(alignment).split())[1:-3])
                    for hsp in alignment.hsps:
                        if not "-" in hsp.query and not "-" in hsp.sbjct:
                            al = Aligner.align(hsp.query, hsp.sbjct)
                            if len(al) > 0:
                                al = al[0]
                                self.alignments[nid][species] = al
                                if al[2] > maxscore[nid]:
                                    self.maximum_alignment[nid] = species
                                    maxscore[nid] = al[2]
        f.close()

    def computeR(self, a=26, k=4.86936):
        '''
        Compute TCR-recognition probabilities for each neoantigen.
        '''
        # iterate over all neoantigens
        for i in self.alignments:
            # energies of all bound states of neoantigen i
            bindingEnergies = list(map(lambda el: -k * (a - el[2]), self.alignments[i].values()))
            # partition function, over all bound states and an unbound state
            lZ = Aligner.logSum(bindingEnergies + [0])
            lGb = Aligner.logSum(bindingEnergies)
            R = exp(lGb - lZ)
            self.Ri[i] = R

    def getR(self, i):
        '''
        Return precomputed R value for a given neoantigen i.
        '''
        if i in self.Ri:
            return self.Ri[i]
        return 0.