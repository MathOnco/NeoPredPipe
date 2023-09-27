#!/usr/bin/env python

'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
Contributions from: Eszter Lakatos

Adapted from Marta Luksza (see RecognitionPotential.md)
'''

class Neoantigen(object):
    '''
    classdocs
    '''
    L = 1.  # default concentration of peptides
    M = 1.  # default concentration of mutant peptides
    W = 1.  # default concentration of wildtype peptides

    WEPS = 0.0003
    WTCAP = float("inf")

    HYDROPHOBIC_RESIDUES = "AILMFWYV"
    WEIRD_RESIDUES = "CGP"

    AGGR = "MAX"
    KDnormalize = 1
    AGGRnum = float("inf")

    @staticmethod
    def residueChangeClass(res1, res2):
        code = ""
        if res1 in Neoantigen.HYDROPHOBIC_RESIDUES:
            code += "H"
        elif res2 in Neoantigen.WEIRD_RESIDUES:
            code += "W"
        else:
            code += "N"
        if res2 in Neoantigen.HYDROPHOBIC_RESIDUES:
            code += "H"
        elif res2 in Neoantigen.WEIRD_RESIDUES:
            code += "W"
        else:
            code += "N"
        return code

    def __init__(self, params, indels):
        '''
        Constructor
        '''
        pparams = params
        if len(params) == 9:
            pparams.append("1")
        [nid, mid, sample, wtPeptide, mtPeptide, allele, wtScore, mtScore, HLA, chopscore] = params[:10]
        self.id = int(nid)
        self.mid = mid
        self.sample = sample
        self.wtPeptide = wtPeptide
        self.mtPeptide = mtPeptide

        if indels==False:
            try:
                [res1, res2] = filter(lambda el: el[0] != el[1], zip(self.wtPeptide, self.mtPeptide))[0]
            except TypeError:
                [res1, res2] = list(filter(lambda el: el[0] != el[1], zip(self.wtPeptide, self.mtPeptide)))[0]

            self.residueChange = Neoantigen.residueChangeClass(res1, res2)
            self.position = list(filter(lambda el: el[1], map(lambda i: [i, self.mtPeptide[i] != self.wtPeptide[i]], range(0, len(self.wtPeptide)))))
            self.position = self.position[0][0] + 1
        else:
            self.residueChange = "-"
            self.position = "FS"
        self.allele = allele
        self.HLA = HLA
        self.chopscore = int(chopscore)
        self.potential = -1e10
        try:
            self.wtkD = min(Neoantigen.WTCAP, float(wtScore))
            self.kD = float(mtScore)
            self.setA()
        except:
            self.kD = Neoantigen.INF
            self.wtkD = Neoantigen.INF
            self.A = 1.
            
        self.expr = None
        if len(params) == 11:
            self.expr = params[10]

    def getSampleName(self):
        return self.sample

    def correctWT(self):
        '''
        Corrects large wildtype kD dissociation constants
        '''
        kd = self.wtkD
        prb = Neoantigen.W / (Neoantigen.W + kd)
        pru = kd / (Neoantigen.W + kd)
        eps = Neoantigen.WEPS
        prb += eps
        pru += eps
        z = 1 + 2 * eps
        prb /= z
        pru /= z
        return pru / prb

    def getWeight(self):
        '''
        Returns 0 for neoantigens that mutated from a non-hydrophobic residues on position 2 or 9;
        these are excluded from analysis. Returns 1 for all other neoantigens
        '''
        w = 1
        if self.residueChange[0] != "H" and (self.position == 2 or self.position == 9):
            w = 0
        return w

    def setA(self):
        '''
        Computes MHC amplitude A
        '''
        self.A = Neoantigen.M / self.kD * self.correctWT()

    def getA(self):
        '''
        Return MHC amplitude A
        '''
        return self.A
