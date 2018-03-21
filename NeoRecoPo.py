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
from NeoClass import Neoantigen
from NeoAlign import Aligner

def Parser():
    # get user variables
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dirty", dest="Dirty", default=True, action='store_false', help="Flag to keep intermediate files. Default: False. Note: Use for debugging.")
    parser.add_argument("-o", '--neoreco_out', dest="neorecoOut", default="", type=str, help="Output Directory.")
    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument("-i", '--neopred_in', dest="neoPredIn", default=None, type=str,help="Input neoantigen predictions. Example: -I ./AllSamples.neoantigens.txt")
    Options = parser.parse_args()  # main user args
    if not Options.neoPredIn:
        parser.error("Some of the required arguments were not provided. Please check required arguments.")
    if len(Options.neorecoOut)!=0:
        if Options.neorecoOut[len(Options.neorecoOut)-1]!='/':
            Options.neorecoOut = Options.neorecoOut+'/'
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

def main():
    # Pull information about usr system files
    localpath = os.path.abspath(__file__).replace('NeoRecoPo.py', '')  # path to scripts working directory
    Config = configparser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    Options = Parser()

    tmpOut = '%s/%s/NeoRecoTMP/'%(Options.neorecoOut, localpath)
    if os.path.isdir(tmpOut)==False:
        os.system('mkdir %s'%(tmpOut))

    if Options.Dirty:
        os.system('rm -r %s'%(tmpOut))


if __name__=="__main__":
    main()