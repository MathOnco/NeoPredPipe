'''
@author: Eszter Lakatos
'''

import unittest
import subprocess
import os

from postprocessing import DefineGenotypeFormat, ProcessPepmatch
from process_expression import BuildGeneIDTable
from hla_preprocess import processHLAminerFile, readInHLA2hlaminer, readInHLA2hlahd, composeHLA2File, ConstructAlleles, ConstructAlleles_typeII

class MyTestCase(unittest.TestCase):
    def test_build_expression_ids(self):
        nmID = "NM_025229"
        geneID = "ENSG00000101251"
        transcriptID = "ENST00000284951"
        uscsID = "uc061vme.1"
        self.assertEqual( (geneID,transcriptID,uscsID), (BuildGeneIDTable('./','ensembl_gene')[nmID], BuildGeneIDTable('./','ensembl_transcript')[nmID], BuildGeneIDTable('./','uscs')[nmID]) )


    def test_filter_read_in_pepmatch(self):
        pmfileName = 'test/Test_pepmatch.out'
        eplines = ['6\tHLA-C*07:02\tTLASKITGM\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline195_NM_0025\t0.1744960\t1.6035',
                   '6\tHLA-C*07:02\tASKITGMLL\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline195_NM_0025\t0.1744960\t1.6035\t<=\tWB',
                   '6\tHLA-C*07:02\tSKITGMLLE\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline195_NM_0025\t0.1744960\t1.6035\t<=\tWB',
                   '6\tHLA-C*07:02\tRLFPLIQAL\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline196_NM_0025\t0.1744960\t1.6035\t<=\tWB']

        appendedlines = ['6\tHLA-C*07:02\tTLASKITGM\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline195_NM_0025\t0.1744960\t1.6035\t1',
                         '6\tHLA-C*07:02\tASKITGMLL\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline195_NM_0025\t0.1744960\t1.6035\t<=\tWB\t0',
                         '6\tHLA-C*07:02\tSKITGMLLE\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline195_NM_0025\t0.1744960\t1.6035\t<=\tWB\t0',
                         '6\tHLA-C*07:02\tRLFPLIQAL\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline196_NM_0025\t0.1744960\t1.6035\t<=\tWB\t1']
        self.assertEqual(appendedlines, ProcessPepmatch(pmfileName, eplines))

    def test_genotypeformat_ad(self):
        line_ad = "line3\tnonsynonymous SNV\tPRAMEF20:NM_001099852:exon2:c.G247A:p.D83N,\tchr1\t13743058\t13743058\tG\tA\t0.1667\t19.4939\t26\tchr1\t13743058\t.\tG\tA\t19.4939\tPASS\tECNT=1;HCNT=22;MAX_ED=.;MIN_ED=.;NLOD=27.62;TLOD=10.35\tGT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:QSS:REF_F1R2:REF_F2R1\t0/0:93,0:0.00:0:0:.:2406,0:49:44\t0/1:57,6:0.081:5:1:0.167:1457,187:32:25"
        self.assertEqual(('alldepths', 1), DefineGenotypeFormat(line_ad))

    def test_genotypeformat_a(self):
        line_a = "line3\tnonsynonymous SNV\tPRAMEF20:NM_001099852:exon2:c.G247A:p.D83N,\tchr1\t13743058\t13743058\tG\tA\t0.1667\t19.4939\t26\tchr1\t13743058\t.\tG\tA\t19.4939\tPASS\tNS=3;DISTR=|G|AG|G|;SB=1.0000	GT:A:GQ:SS:BCOUNT:DP\t0/0:G:100.0000:0:0,0,18,0:18\t0/1:AG:19.4939:2:4,0,15,0:19\t0/0:G:100.0000:0:0,0,26,0:26"
        self.assertEqual(('allele', 1), DefineGenotypeFormat(line_a))

    def test_genotypeformat_nv(self):
        line_nv = "line3\tnonsynonymous SNV\tPRAMEF20:NM_001099852:exon2:c.G247A:p.D83N,\tchr1\t13743058\t13743058\tG\tA\t0.1667\t19.4939\t26\tchr1\t13743058\t.\tG\tA\t19.4939\tMQ;badReads\tAC=5;AF=0.500;AN=10;BRF=0.97;FR=0.4556;HP=6;HapScore=2;MGOF=39;MMLQ=26;MQ=0.36;NF=4;NR=2;PP=111;QD=24.2952020912;SC=CAGATAGTGGAGGGGCTTACA;SbPval=0.65;Source=Platypus;TC=14;TCF=10;TCR=4;TR=6;WE=621651;WS=621636;set=FilteredInAll\tGT:GOF:GQ:NR:NV:PL\t0/1:4:15:2:1:33,0,15"
        self.assertEqual(('numvarreads', 4), DefineGenotypeFormat(line_nv))

    def test_genotypeformat_freq(self):
        line_freq = "line3\tnonsynonymous SNV\tPRAMEF20:NM_001099852:exon2:c.G247A:p.D83N,\tchr1\t13743058\t13743058\tG\tA\t0.1667\t19.4939\t26\tchr1\t13743058\t.\tG\tA\t19.4939\tPASS\tECNT=1;HCNT=22;MAX_ED=.;MIN_ED=.;NLOD=27.62;TLOD=10.35\tGT:GQ:DP:RD:AD:FREQ:DP4\t0/0:.:21:21:0:0%:20,1,0,0\t0/1:.:28:23:5:17.86%:22,1,5,0"
        self.assertEqual(('varscanfreq',5), DefineGenotypeFormat(line_freq))

    def test_genotypeformat_gt(self):
        line_gt = "line3\tnonsynonymous SNV\tPRAMEF20:NM_001099852:exon2:c.G247A:p.D83N,\tchr1\t13743058\t13743058\tG\tA\t0.1667\t19.4939\t26\tchr1\t13743058\t.\tG\tA\t19.4939\tPASS\tECNT=1;HCNT=22;MAX_ED=.;MIN_ED=.;NLOD=27.62;TLOD=10.35\tGT:IGT:DP:DP4:BCOUNT:GQ:JGQ:VAQ:BQ:MQ:AMQ:SS:SSC\t0/0:0/0:8:5,3,0,0:0,0,0,8:51:19:0:26:12:12:0:.\t0/1:0/1:7:2,3,1,1:2,0,0,5:12:19:12:31,28:20:30,15:2:19"
        self.assertEqual(('genotype',0), DefineGenotypeFormat(line_gt))

    def test_genotypeformat_strelka(self):
        line_strelka = "line984\tsynonymous SNV\tRAB25:NM_020387:exon2:c.C114T:p.S38S,\tchr1\t156065981\t156065981\tC\tT\t.\t.\t.\tchr1\t156065981\t.\tC\tT\t.\tPASS\tSOMATIC;QSS=56;TQSS=1;NT=ref;QSS_NT=56;TQSS_NT=1;SGT=CC->CT;DP=267;MQ=60.00;MQ0=0;ReadPosRankSum=-1.09;SNVSB=0.00;SomaticEVS=10.10\tDP:FDP:SDP:SUBDP:AU:CU:GU:TU\t162:0:0:0:0,0:162,163:0,0:0,0\t103:0:0:0:0,0:99,100:0,0:4,4"
        self.assertEqual(('strelka',5), (DefineGenotypeFormat(line_strelka)[0],DefineGenotypeFormat(line_strelka)[1]['CU']))

    def test_hla_format_typeI(self):
        hlas = ['hla_a_03_01_27', 'hla_a_01_01_01_01', 'hla_b_07_02_09', 'hla_a_33_03_03q', 'hla_a_03_01_01_02n','hla_b_39_01_01_02l','hla_b_82_02', 'NA']
        FilePath = '.'
        patID = 'hla_test'

        self.assertEqual( set(['HLA-A33:03','HLA-B82:02','HLA-A01:01','HLA-A03:01','HLA-B39:01','HLA-B07:02']), set(ConstructAlleles(hlas,FilePath,patID) ))

    def test_hla_format_typeII(self):
        hlas = ['DRB1*12:02P', 'DPA1*01:03P','DPB1*90:01P','DRB1*01:01:01G', 'DPA1*02:04','NA', 'DQA1*05:01P', 'DQB1*02:02:03:01', 'DQB1*03:39']
        FilePath = '.'
        patID = 'hla_test'

        self.assertEqual( set(['HLA-DQA10501-DQB10202','HLA-DPA10103-DPB19001','DRB1_0101','DRB1_1202','HLA-DPA10204-DPB19001']), set(ConstructAlleles_typeII(hlas,FilePath,patID)))

    def test_hla_process_hlaminer_typeII(self):
        if os.path.isfile("./test/hla-II/Test_hlaminer/HLAminer_processed.txt"):
            os.system("rm ./test/hla-II/Test_hlaminer/HLAminer_processed.txt")

        hlas1 = ['DPA1*01:03P','DPA1*03:01P','DPB1*02:01P','DPB1*04:01P','DQA1*01:02P','DQA1*05:01P','DQB1*06:02P','DQB1*03:01P','DRB1*09:01P','DRB1*07:01P']
        hlas2 = ['DRB1*13:01:01','DRB1*03:01:01','DQA1*05:01:01','DQA1*01:03:01','DQB1*06:03:01','DQB1*02:01:01','DPA1*01:03:01','DPB1*04:02:01','DPB1*04:01:01']

        hlaDict = composeHLA2File("./test/hla-II")

        self.assertEqual( (hlas1,hlas2), (hlaDict['Test_hlaminer'],hlaDict['Test_hlahd']))

    def test_main_multiallele(self):
        if os.path.isfile("./test/Test_multiAllele.neoantigens.unfiltered.txt"):
            os.system("rm ./test/Test_multiAllele.*")

        cmd = ['python', 'NeoPredPipe.py', '-I', './test/vcfs/', '-H', './test/hlatypes_multiallele.txt', '-o', './test/',
               '-n', 'Test_multiAllele', '-c', '0', '1', '2', '3', '-E', '9', '-a']

        runcmd = subprocess.Popen(cmd)
        runcmd.wait()
        with open('test/Test_multiAllele.neoantigens.unfiltered.txt', 'r') as testof:
            oflines = testof.readlines()
        self.assertEqual( (['1', '1', '0', '0'],['1','1','1','0']) , (oflines[0].split('\t')[1:5], oflines[9].split('\t')[1:5]))

    def test_main_multiple(self):
        if os.path.isfile("./test/Test_platypus.neoantigens.txt"):
            os.system("rm ./test/Test_platypus.*")

        cmd = ['python', 'NeoPredPipe.py', '-I', './test/vcfs/', '-H', './test/hlatypes.txt', '-o', './test/',
               '-n', 'Test_platypus', '-c', '0', '1', '2', '3', '4', '-E', '8', '-d', '-m', '-x', './test/expression.txt' ]

        runcmd = subprocess.Popen(cmd)
        runcmd.wait()
        with open('test/Test_platypus.neoantigens.txt', 'r') as testof:
            oflines = testof.readlines()
        self.assertEqual( ['1', '1', '0', '1', '0'] , oflines[0].split('\t')[1:6])
        
    def test_main_multiple_summaries(self):
        with open('test/Test_platypus.neoantigens.summarytable.txt', 'r') as testsum:
            sumlines = testsum.readlines()
        summary = sumlines[1].rstrip('\n').split('\t')
        #self.assertEqual( (['3','3','2','2','2'], ['1','0','0','0','1','1']), (summary[4:9], summary[22:])) #true for EL
        self.assertEqual( (['3','3','2','1','2'], ['0','0','0','0','2','1']), (summary[4:9], summary[22:]))

    def test_main_peptide_checking(self):
        with open('test/Test_platypus.neoantigens.txt', 'r') as testof:
            oflines = testof.readlines()
       # self.assertEqual( ('0', '1'), (oflines[1].rstrip('\n').split('\t')[-1], oflines[2].rstrip('\n').split('\t')[-1])) #true for EL
        self.assertEqual( ('1', '1'), (oflines[1].rstrip('\n').split('\t')[-1], oflines[2].rstrip('\n').split('\t')[-1]))

    def test_main_read_expression(self):
        with open('test/Test_platypus.neoantigens.txt', 'r') as testof:
            oflines = testof.readlines()
        self.assertEqual( ('NA', '0.12'), (oflines[0].rstrip('\n').split('\t')[-18], oflines[1].rstrip('\n').split('\t')[-18]))


    def test_main_recopo(self):
        if os.path.isfile("./test/PredictedRecognitionPotentials.txt"):
            os.system("rm ./test/PredictedRecognitionPotentials.txt")
        cmd = ['python', 'NeoRecoPo.py', '-i', './test/Test_platypus.neoantigens.txt', '-f', './test/fastaFiles/', '-o', './test/']
        runcmd = subprocess.Popen(cmd)
        runcmd.wait()
        with open('test/PredictedRecognitionPotentials.txt', 'r') as testof:
            oflines = testof.readlines()
        self.assertEqual(['1', 'line3_NM_001005', 'Test_platypus', '3', 'HN', 'KPRHYLTI', 'KPLHYLTI', 'B0702', '0.12', '7.54006501848'], oflines[1].split('\t')[:-3])

    def test_main_single_region(self):
        if os.path.isfile("test/Test_single.neoantigens.summarytable.txt"):
            os.system("rm ./test/Test_single.*")
        cmd = ['python', 'NeoPredPipe.py', '-I', './test/vcfs/', '-H', './test/hlatypes.txt', '-o', './test/',
               '-n', 'Test_single', '-E', '8' ]
        runcmd = subprocess.Popen(cmd)
        runcmd.wait()
        with open('test/Test_single.neoantigens.summarytable.txt', 'r') as testof:
            oflines = testof.readlines()
        self.assertEqual( ['3', '2', '1'] , oflines[1].rstrip('\n').split('\t')[1:])

    def test_main_strelka(self):
        if os.path.isfile("test/Test_strelka.neoantigens.txt"):
            os.system("rm ./test/Test_strelka.*")
        cmd = ['python', 'NeoPredPipe.py', '-I', './test/vcfs/', '-H', './test/hlatypes_strelka.txt', '-o', './test/',
               '-n', 'Test_strelka', '-E', '8' , '-c', '1','2','3','--manualproc']
        runcmd = subprocess.Popen(cmd)
        runcmd.wait()
        with open('test/Test_strelka.neoantigens.txt', 'r') as testof:
            oflines = testof.readlines()
        self.assertEqual( (['0','1','1'],['1','1','1']) , (oflines[0].split('\t')[1:4],oflines[1].split('\t')[1:4]) )


if __name__ == '__main__':
    import sys
    print(sys.version)
    unittest.main()
