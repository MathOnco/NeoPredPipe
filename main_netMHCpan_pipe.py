#!/usr/bin/env python2.7

import sys
import os
import argparse
import ConfigParser

# get user variables
parser = argparse.ArgumentParser()
parser.add_argument("InputCSV", type=str, help=".csv file with Input VCF, normal BAM, semi-colon separated epitope lengths, Polysolver HLA Type File")
##parser.add_argument("I", type=str, help="Input vcf files location")
##parser.add_argument("E", type=str, help="Comma separated epitope lengths for prediction (e.g. 8,9,10)")
##parser.add_argument("H", type=str, help="Polysolver HLA Type file")
parser.add_argument("OutputDir", type=str, help="Output Directory Path")
args = parser.parse_args() # main user args

def ConfigSectionMap(section):
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

    def __init__(self):
        self.Number = None
        self.test = self.setit()

    def setit(self):
        self.Number=100

def PrepClasses():
    t = []
    for i in range(0,1):
        t.append(Sample)
    return t

if __name__=="__main__":
    # Pull information about usr system files
    localpath = os.path.abspath(__file__).rstrip('main_netMHCpan_pipe.py')  # path to scripts working directory
    Config = ConfigParser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    annPaths =ConfigSectionMap(Config.sections()[0]) # get annovar script paths

    t = PrepClasses()
    print t[0].test

