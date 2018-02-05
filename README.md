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
1. You can see the options associated by running the following:
```bash
python ./main_netMHCpan_pipe.py
```