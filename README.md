## Summary

This tool allows a user to process neoantigens predicted from vcf files using ANNOVAR and netMHCpan.

## Dependencies
##### Note: Should be compatible on Darwin and Linux systems, not Windows.

1. Python == 2.7 (Built using Python 2.7.13, not compatible with python 3 due to OS processes)
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

## Getting Started
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
                              [-I VCFDIR] [-H HLAFILE] [-o OUTPUTDIR] [-pp]
                              [-m] [-c COLREGIONS [COLREGIONS ...]] [-a] [-t]

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

Required arguments:
  -I VCFDIR             Input vcf file directory location. Example: -I
                        ./Example/input_vcfs/
  -H HLAFILE            HLA file for vcf patient samples.
  -o OUTPUTDIR          Output Directory Path

Post Processing Options:
  -pp                   Flag to perform post processing. Default=True.
  -m                    Specifies if the vcf is a multiregion sample. Default:
                        False.
  -c COLREGIONS [COLREGIONS ...]
                        Columns of regions within vcf that are not normal
                        within a multiregion vcf file. 0 is normal in test
                        samples. Can handle different number of regions per
                        vcf file.
  -a                    Flag to not filter neoantigen predictions and keep all
                        regardless of prediction value.
  -t                    Flag to turn off summary table.
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

## Example of how to run the script using the provided example vcf files and corresponding hla types
```bash
python main_netMHCpan_pipe.py -I ./Example/input_vcfs -H ./Example/HLAtypes/hlatypes.txt -o ./ -m -c 1 2 -E 8 9 10
```

## Data post processing
1. Post processing is turned on by default. If you want it turned off set the '-pp' flag.
2. The output files will yield files with the following information:
   - A file containing the neoantigen predictions with appropriate identifer information
   - TODO add an example of output with what is in each column
   - A file containing summaries of the neoantigen burdens in each sample (and regions if multiregion).
3. If there are not multiple regions from a single patient the resulting summary table will appear as follows:

| Patient | Total | Total_WB | Total_SB |
|  --- |  --- |  --- |  --- |

4. If multiple regions are specified then the output will look as follows:

| Patient | Total | Total_WB | Total_SB | Total_Region_1 | Total_Region_n | Total_WB_Region_1 | Total_WB_Region_n | Total_SB_Region_1 | Total_SB_Region_n | Clonal | Subclonal | Shared | Clonal_WB | Clonal_SB | Subclonal_WB | Subclonal_SB | Shared_WB | Shared_SB |
|  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |
