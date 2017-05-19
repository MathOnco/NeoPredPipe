#!/usr/bin/env python2.7

import sys
import os
import argparse
import ConfigParser

# get user variables
parser = argparse.ArgumentParser()
parser.add_argument("I", type=str, help="Input vcf files location")
parser.add_argument("E", type=str, help="Comma separated epitope lengths for prediction (e.g. 8,9,10)")
parser.add_argument("H", type=str, help="Polysolver HLA Type file")
parser.add_argument("O", type=str, help="Output Directory Path")
args = parser.parse_args() # main user args
localpath = os.path.abspath(__file__).rstrip('main_netMHCpan_pipe.py') # path to scripts working directory


class MyParser(ConfigParser.ConfigParser):

    def as_dict(self):
        d = dict(self._sections)
        for k in d:
            d[k] = dict(self._defaults, **d[k])
            d[k].pop('__name__', None)
        return d

MyParser.as_dict()