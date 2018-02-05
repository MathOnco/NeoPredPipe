## Getting Started

This tool allows a user to process neoantigens predicted from vcf files using ANNOVAR and netMHCpan.

## Dependencies

1. Python >= 2.7 (Built using Python 2.7.13)
   - argparse
   - ConfigParser
2. ANNOVAR
3. hg19 reference genome
4. netMHCpan

## Running
1. Clone the repository:
```bash
git clone https://github.com/rschenck/NeoPredPipe.git
```
2. You can see the options associated by running the following:
```bash
python ./main_netMHCpan_pipe.py --help
```
   - Which produces the following:
```bash
usage: main_netMHCpan_pipe.py [-h] [-E EPITOPES [EPITOPES ...]] [-I VCFDIR]
                              [-H HLAFILE] [-o OUTPUTDIR]

optional arguments:
  -h, --help            show this help message and exit
  -E EPITOPES [EPITOPES ...], --epitopes EPITOPES [EPITOPES ...]
                        Epitope lengths for predictions. Default: 8 9 10

Required arguments:
  -I VCFDIR             Input vcf file directory location. Example: --VCFDir
                        ./Example/input_vcfs/
  -H HLAFILE            HLA file for vcf patient samples.
  -o OUTPUTDIR          Output Directory Path
```

## Input files
1. VCF file. A standard vcf file with a patient identifier as the title of the .vcf.
2. An hla file with the following format:
   - Note, patient identifier in the rows should match that in the vcf.
   - Headers are not required but the data should match the format in the table.
Patient | HLA-A_1 | HLA-A_2 | HLA-B_1 | HLA-B_2 | HLA-C_1 | HLA-C_2 |
 --- |  --- |  --- |  --- |  --- |  --- |  ---  |
12-N | hla_a_31_01_02 | hla_a_02_01_80 | hla_b_40_01_02 | hla_b_50_01_01 | hla_c_03_04_20 | hla_c_06_02_01_02 |
XN | hla_a_01_01_01_01 | NA | hla_b_07_02_01 | NA | hla_c_01_02_01 | NA |


## UNDER DEVELOPMENT