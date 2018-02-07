## Summary

This tool allows a user to process neoantigens predicted from vcf files using ANNOVAR and netMHCpan.

## Dependencies
##### Note: Should be compatible on Darwin and Linux systems, not Windows.

1. Python == 2.7 (Built using Python 2.7.13, not compatible with python 3 yet)
   - biopython == 1.70
2. ANNOVAR
   - Can be downloaded [here](http://annovar.openbioinformatics.org/en/latest/user-guide/download/).
   - ANNOVAR hg19_refGene
   - ANNOVAR hg19_refGeneMrna
   - Other reference builds can be used. Simply change the usr_path.ini file to the appropriate reference (see below).
     - Make sure to use the same one used to call variants.
4. netMHCpan
   - Using [netMHCpan-4.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan) for all tests of this pipeline.
   - Follow their steps for installation on your platform.

## Installing and preparing environment
1. Clone the repository:
```bash
git clone https://github.com/rschenck/NeoPredPipe.git
```
2. Configure the 'usr_path.ini' file for your environment.
   - All paths within the annovar header should be where you installed annovar.
   - Only one path is needed to the netMHCpan executible under netMHCpan
##### Note: You need to provide the absolute path.
3. You can see the options associated by running the following:
```bash
python ./main_netMHCpan_pipe.py --help
```
   - Which produces the following:
```bash
usage: main_netMHCpan_pipe.py [-h] [-E EPITOPES [EPITOPES ...]] [-l] [-d] [-r]
                              [-p] [-I VCFDIR] [-H HLAFILE] [-o OUTPUTDIR]
                              [-n OUTNAME] [-pp]
                              [-c COLREGIONS [COLREGIONS ...]] [-a] [-t]

optional arguments:
  -h, --help            show this help message and exit
  -E EPITOPES [EPITOPES ...], --epitopes EPITOPES [EPITOPES ...]
                        Epitope lengths for predictions. Default: 8 9 10
  -l                    Specifies whether to delete the ANNOVAR log file.
                        Default: True. Note: Use for debugging.
  -d                    Specified whether to delete intermediate files created
                        by program. Default: True. Note: Set flag to resume
                        job.
  -r, --cleanrun        Specify this alone with no other options to clean-up a
                        run. Be careful that you mean to do this!!
  -p, --preponly        Prep files only without running neoantigen
                        predictions. The prediction step takes the most time.

Required arguments:
  -I VCFDIR             Input vcf file directory location. Example: -I
                        ./Example/input_vcfs/
  -H HLAFILE            HLA file for vcf patient samples.
  -o OUTPUTDIR          Output Directory Path
  -n OUTNAME            Name of the output file for neoantigen predictions

Post Processing Options:
  -pp                   Flag to perform post processing. Default=True.
  -c COLREGIONS [COLREGIONS ...]
                        Columns of regions within vcf that are not normal
                        within a multiregion vcf file after the format field.
                        Example: 0 is normal in test samples, tumor are the
                        other columns. Program can handle different number of
                        regions per vcf file.
  -a                    Flag to not filter neoantigen predictions and keep all
                        regardless of prediction value.
  -t                    Flag to turn off a neoantigen burden summary table.
                        Default=True.
```

## Input files
1. VCF file. A standard vcf file with a patient identifier as the title of the .vcf.
2. An hla file with the following tab delimited format:
   - Note, patient identifier in the rows must match that preceding *.vcf
   - Headers are not required but the data should match the format in the table.
   - 'NA' is used when the HLA typing predicts the same HLA subtype for A, B, or C.
   - The program will search for the appropriate allele within netMHCpan alleles list, but care should be taken to ensure accuracy.

| Patient | HLA-A_1 | HLA-A_2 | HLA-B_1 | HLA-B_2 | HLA-C_1 | HLA-C_2 |
|  --- |  --- |  --- |  --- |  --- |  --- |  ---  |
| test1 | hla_a_31_01_02 | hla_a_02_01_80 | hla_b_40_01_02 | hla_b_50_01_01 | hla_c_03_04_20 | hla_c_06_02_01_02 |
| test2 | hla_a_01_01_01_01 | NA | hla_b_07_02_01 | NA | hla_c_01_02_01 | NA |

## Run Using Example .vcf files
```bash
# Run the Pipeline to only prepare the input files. Can be best to run this independent if working on a cluster.
python main_netMHCpan_pipe.py --preponly -I ./Example/input_vcfs -H ./Example/HLAtypes/hlatypes.txt -o ./ -n TestRun -c 1 2 -E 8 9 10

# Run the Pipeline
python main_netMHCpan_pipe.py -I ./Example/input_vcfs -H ./Example/HLAtypes/hlatypes.txt -o ./ -n TestRun -c 1 2 -E 8 9 10
```

## Data post processing
1. Post processing is turned on by default. If you want it turned off set the '-pp' flag.
2. The output files will yield files with the following information:
   - A file containing the neoantigen predictions with appropriate identifier information and heterogeneity if multiregion.
   - A file containing summaries of the neoantigen burdens in each sample (and regions if multiregion).

## Output Format
1. The primary output file of neoantigens has the following format:
   - Stuff

| Sample |  R1 |  R2 |  R3 |  Line |  chr |  allelepos |  ref |  alt |  GeneName:RefSeqID |  pos |  hla |  peptide |  core |  Of |  Gp |  Gl |  Ip |  Il |  Icore |  Identity |  Score |  Aff |  Rank |  Candidate | BindLevel |
| --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |
| test1 | 0 | 1 | 0 | line16 | chr1 | 153914523 | G | C | DENND4B:NM_014856 | 3 | HLA-B*40:01 | SERQAGAL | SERQAG-AL | 0 | 0 | 0 | 6 | 1 | SERQAGAL | line16_NM_01485 | 0.33670 | 1308.7 | 1.30 | <= | WB |
| test1 | 1 | 1 | 0 | line8 | chr1 | 53608000 | C | T | SLC1A7:NM_001287597,SLC1A7:NM_001287595,SLC1A7:NM_006671,SLC1A7:NM_001287596 | 2 | HLA-C*06:02 | LGFFLRTRHL | LFFLRTRHL | 0 | 1 | 1 | 0 | 0 | LGFFLRTRHL | line8_NM_001287 | 0.24655 | 3470.8 | 1.20 | <= | WB |
| test2 | 1 | 0 | 0 | line34 | chr1 | 248402593 | C | A | OR2M4:NM_017504 | 6 | HLA-C*01:02 | VMAYERYVAI | VAYERYVAI | 0 | 1 | 1 | 0 | 0 | VMAYERYVAI | line34_NM_01750 | 0.14917 | 9954.7 | 1.50 | <= | WB |
| test2 | 1 | 1 | 0 | line51 | chr2 | 240982213 | C | G | PRR21:NM_001080835 | 2 | HLA-C*01:02 | FTHGPSSTPL | FTHPSSTPL | 0 | 3 | 1 | 0 | 0 | FTHGPSSTPL | line51_NM_00108 | 0.22570 | 4349.1 | 0.40 | <= | SB |
| test2 | 1 | 1 | 0 | line51 | chr2 | 240982213 | C | G | PRR21:NM_001080835 | 7 | HLA-C*01:02 | SSTPLHPCPF | STPLHPCPF | 0 | 1 | 1 | 0 | 0 | SSTPLHPCPF | line51_NM_00108 | 0.13137 | 12068.7 | 2.00 | <= | WB |

2. If there are not multiple regions from a single patient the resulting summary table will appear as follows:

| Patient | Total | Total_WB | Total_SB |
|  --- |  --- |  --- |  --- |

4. If multiple regions are specified then the output will look as follows (scroll left or right to view all):

| Patient | Total | Total_WB | Total_SB | Total_Region_1 | Total_Region_n | Total_WB_Region_1 | Total_WB_Region_n | Total_SB_Region_1 | Total_SB_Region_n | Clonal | Subclonal | Shared | Clonal_WB | Clonal_SB | Subclonal_WB | Subclonal_SB | Shared_WB | Shared_SB |
|  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |
