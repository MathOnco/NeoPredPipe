## Summary

This tool allows a user to process neoantigens predicted from single- or multi-region vcf files using ANNOVAR and netMHCpan.

Once the primary pipeline is ran, the user is then able to perform Neoantigen recognition potential as described by [Marta Luksza et al., Nature 2017](https://www.nature.com/articles/nature24473).
To perform the neoantigen recognition potential please [click here](RecognitionPotential.md). Finally, you can visualize these predictions using [NeoPredViz](NeoPredViz.md).

## Citation

NeoPredPipe has now been published in BMC Bioinformatics. Please cite:
- Schenck, Ryan O. and Lakatos, Eszter and Gatenbee, Chandler and Graham, Trevor A. and Anderson, Alexander R. A. [NeoPredPipe: high-throughput neoantigen prediction and recognition potential pipeline](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2876-4). BMC Bioinformatics. 2019. 20:264.

## Dependencies
##### Note: Should be compatible on Darwin and Linux systems, not Windows.

1. Python == 2.7 (Built using Python 2.7.13, not compatible with python 3 yet)
   - biopython == 1.70
2. ANNOVAR
   - Can be downloaded [here](http://annovar.openbioinformatics.org/en/latest/user-guide/download/).
   - ANNOVAR reference files, for example: hg19_refGene.txt and hg19_refGeneMrna.fa
   - Other reference files/builds can be used. Simply change the usr_path.ini file to the appropriate reference (see below).
     - Make sure to use the same one used to call variants.
   - **NOTE: For indel predictions, we highly recommend to use the 2018-04-16 or later release of ANNOVAR, as earlier versions do not provide the appropriate support for protein-elongating frameshift mutations.**
4. netMHCpan or netMHCIIpan
   - The pipeline has been primarily developed to support [netMHCpan-4.0](http://www.cbs.dtu.dk/services/NetMHCpan-4.0/). It is also compatible with netMHCpan-4.1, but have not been tested with this latest version.
   - Support for MHC-II prediction using [netMHCIIpan-3.2](http://www.cbs.dtu.dk/services/NetMHCIIpan-3.2/) is also available, but note that this functionality has not been thoroughly tested.  
   - To install either software, choose your platform and follow instructions [here](https://services.healthtech.dtu.dk/software.php).
5. PeptideMatch (Only necessary if one wishes to check predicted epitopes for novelty against a reference proteome.)
   - Requires Java.
   - The runnable jar is available [here](https://research.bioinformatics.udel.edu/peptidematch/commandlinetool.jsp).
   - Download a reference protein sequence in fasta format (e.g. from [Ensembl](ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/pep/) or [UniProt](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/)) and index it according to the Tutorial.
   - **We advise the use of PeptideMatch for indel predictions, to filter out non-frameshift peptides and peptides that are novel to the genomic location, but coincidentally exist elsewhere.**

## Installing and preparing environment
1. Clone the repository:
```bash
git clone https://github.com/MathOnco/NeoPredPipe.git
```
2. Configure the 'usr_path.ini' file for your environment.
   - All paths within the annovar header should be where you installed annovar.
   - Only one path is needed to the netMHCpan executible under netMHCpan
      - If you wish to use netMHCIIpan instead, simply provide that path in this section and the pipeline will process files accordingly
   - If you wish to use PeptideMatch, provide paths for both jar and reference index.
   #### Note: You need to provide the absolute paths for all the above.
3. You can see the options associated by running the following:
```bash
python ./NeoPredPipe.py --help
```
   - Which produces the following:
```bash
usage: NeoPredPipe.py [-h] [-E EPITOPES [EPITOPES ...]] [-l] [-d] [-r] [-p]
                      [--EL] [-I VCFDIR] [-H HLAFILE] [-o OUTPUTDIR]
                      [-n OUTNAME] [-pp] [-c COLREGIONS [COLREGIONS ...]] [-a]
                      [-m] [-x EXPRESSION] [--expmulti] [-t]

optional arguments:
  -h, --help            show this help message and exit
  -E EPITOPES [EPITOPES ...], --epitopes EPITOPES [EPITOPES ...]
                        Epitope lengths for predictions. Default: 8 9 10
  -l                    Specifies whether to delete the ANNOVAR log file.
                        Default: True. Note: Use for debugging.
  -d                    Specifies whether to delete intermediate files created
                        by program. Default: True. Note: Set flag to resume
                        job.
  -r, --cleanrun        Specify this alone with no other options to clean-up a
                        run. Be careful that you mean to do this!!
  -p, --preponly        Prep files only without running neoantigen
                        predictions. The prediction step takes the most time.
  --EL                  Flag to perform netMHCpan predictions with Eluted
                        Ligand option (without the -BA flag). Please note that
                        the output will NOT be compatible with downstream
                        Recognition Potential analysis. Default=False (BA
                        predictions)

Required arguments:
  -I VCFDIR             Input vcf file directory location. Example: -I
                        ./Example/input_vcfs/
  -H HLAFILE            HLA file for vcf patient samples OR directory with
                        patient-specific directories from running POLYSOLVER
                        (see Readme).
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
  -m                    Specifies whether to perform check if predicted
                        epitopes match any normal peptide. If set to True,
                        output is added as a column to neoantigens file.
                        Requires PeptideMatch specified in usr_paths.ini.
                        Default=False
  -x EXPRESSION, --expression EXPRESSION
                        RNAseq expression quantification file(s), if
                        specified, expression information is added to output
                        tables.
  --expmulti            Flag to specify if expression file(s) has information
                        on multiple regions in multiple columns.
                        Default=False.
  -t                    Flag to turn off a neoantigen burden summary table.
                        Default=True.
```

## Input files
1. VCF file(s). A standard vcf file with a patient identifier as the title of the .vcf. Several vcf files can be specified in the same directory.
2. An hla file with the following tab delimited format:
   - Note, patient identifier in the rows must match that preceding *.vcf
   - Patient identifier and HLA types should be separated by tabulators.
   - Headers are not required but the data should match the format in the table.
   - 'NA' is used when the HLA typing predicts the same HLA subtype for A, B, or C.
   - The program will search for the appropriate allele within netMHCpan alleles list, but care should be taken to ensure accuracy.

| Patient | HLA-A_1 | HLA-A_2 | HLA-B_1 | HLA-B_2 | HLA-C_1 | HLA-C_2 |
|  --- |  --- |  --- |  --- |  --- |  --- |  ---  |
| test1 | hla_a_31_01_02 | hla_a_02_01_80 | hla_b_40_01_02 | hla_b_50_01_01 | hla_c_03_04_20 | hla_c_06_02_01_02 |
| test2 | hla_a_01_01_01_01 | NA | hla_b_07_02_01 | NA | hla_c_01_02_01 | NA |

2. (Option 2) In case patient-specific HLA types are not known, but WGS/WES bam files are available for the patient, we recommend the use of _in silico_ haplotyping tool, [POLYSOLVER](https://software.broadinstitute.org/cancer/cga/polysolver). POLYSOLVER can be obtained in a cross-platform, easy-to-use Docker container format, as described [here](https://software.broadinstitute.org/cancer/cga/polysolver_download).
   - Follow the steps on the link above for setup and testing POLYSOLVER-Docker.
   - Modify and execute the code below to run polysolver on a sample with name <patientIdentifer>.
   - After executing this for each sample, simply supply the *absolute* path to the directory "HLAtypes" as the -H input of NeoPredPipe, and your HLA predictions will be automatically read into the program.
   
```bash
# Specify your home folder that has a folder of bamfiles and a folder where you want to save the HLAtyping outputs too
mydir="/path/to/your/home/directory"

# Create a folder for the output of haplotyping
mkdir $mydir/HLAtypes
mkdir $mydir/HLAtypes/<patientIdentifier>

# Change paths to fit the names of your own folders/files (and change hg19 to hg38 if needed) and run the following command to generate outputs in a folder named <patientIdentifier>
docker run -d -P --name $NAME -v $mydir:/data sachet/polysolver:v4 bash /home/polysolver/scripts/shell_call_hla_type /data/<bamfiles>/<patientIdentifier>.bam Unknown 1 hg19 STDFQ 0 /data/HLAtypes/<patientIdentifier>
```

3. (Optional) Expression file(s), specified after the -x flag:
   - Either the path to a single file to be used for all samples (for example values compiled from a reference cohort); or a path to a directory, containing files for each sample, named as <patientIdentifier>.tsv. NeoPredPipe will automatically search for appropriate *.tsv files if a directory is specified.
   - Each file should be tab-delimited, with the first column containing gene identifiers, and the second column containing expression values. If the --expmulti flag is specified, expression files are assumed to hold information on multiple regions, and additional columns are read in as expression data.
   - Direct outputs of RNAseq quantification software, such as HTseqCount and kallisto are supported.
   - Currently supported gene identifier formats: Ensembl gene ID, Ensembl transcript ID, RefSeq transcript ID, UCSC transcript ID.

## Run Using Example .vcf files
```bash
# Run the Pipeline to only prepare the input files. Can be best to run this independently if working on a cluster.
python NeoPredPipe.py --preponly -I ./Example/input_vcfs -H ./Example/HLAtypes/hlatypes.txt -o ./ -n TestRun -c 1 2 -E 8 9 10

# Run the Pipeline
python NeoPredPipe.py -I ./Example/input_vcfs -H ./Example/HLAtypes/hlatypes.txt -o ./ -n TestRun -c 1 2 -E 8 9 10
```

## Data post processing
1. Post processing is turned on by default. If you want it turned off set the '-pp' flag.
2. The output files will yield files with the following information:
   - A file containing the neoantigen predictions with appropriate identifier information and heterogeneity if multiregion.
   - A file containing summaries of the neoantigen burdens in each sample (and regions if multiregion).

## Output Format
1. The primary output file of neoantigens has the following format, separated by tabulators (columns 12-26 are taken from [here](http://www.cbs.dtu.dk/services/NetMHCpan/output.php)):
   - **Sample**: vcf filename/patient identifier
   - **R1**: Region 1 of a multiregion sample, binary for presence (1) or absence (0), regions above the number of regions in the sample (for varying number of biopsies) are indicated by -1. Can be *n* numbers of regions. _Only present in multiregion samples_.
   - **R2**: Region 2 of a multiregion sample, binary for presence (1) or absence (0), regions above the number of regions in the sample (for varying number of biopsies) are indicated by -1. Can be *n* numbers of regions. _Only present in multiregion samples_.
   - **R3**: Region 3 of a multiregion sample, binary for presence (1) or absence (0), regions above the number of regions in the sample (for varying number of biopsies) are indicated by -1. Can be *n* numbers of regions. _Only present in multiregion samples_.
   - **Line**: Line number from the *.avready file (same as the vcf) to identify mutation yielding corresponding neoantigen.
   - **chr**: Chromosome of mutation
   - **allelepos**: Position of the mutation
   - **ref**: Reference base at the position
   - **alt**: Alternative base at the location
   - **GeneName:RefID**: Gene name and RefSeq ID separated by a colon. Multiple genes/RefSeq IDs separated by a comma.
   - **Expression**: Expression value of the gene. Expression values for multiple regions (_if using the -expmulti flag_) are comma-separated. NA for genes that are not found in the corresponding expression file, or for samples without expression information. _Only present if the -x flag is used_.
   - **pos**: Residue number (starting from 0)
   - **hla**: Molecule/allele name
   - **peptide**: Amino acid sequence of the potential ligand
   - **core**: The minimal 9 amino acid binding core directly in contact with the MHC
   - **Of**: The starting position of the Core within the Peptide (if > 0, the method predicts a N-terminal protrusion)
   - **Gp**: Position of the deletion, if any.
   - **Gl**: Length of the deletion.
   - **Ip**: Position of the insertions, if any.
   - **Il**: Length of the insertion.
   - **Icore**: Interaction core. This is the sequence of the binding core including eventual insertions of deletions.
   - **Identity**: Protein identifier, i.e. the name of the Fasta entry.
   - **Score**: The raw prediction score
   - **Binding Affinity**: Predicted binding affinity in nanoMolar units.
   - **Rank**: Rank of the predicted affinity compared to a set of random natural peptides. This measure is not affected by inherent bias of certain molecules towards higher or lower mean predicted affinities. Strong binders are defined as having %rank<0.5, and weak binders with %rank<2. We advise to select candidate binders based on %Rank rather than nM Affinity
   - **Candidate**: Symbol (<=) used to denote a Strong or Week Binder in BindLevel
   - **BindLevel**: (SB: strong binder, WB: weak binder). The peptide will be identified as a strong binder if the % Rank is below the specified threshold for the strong binders, by default 0.5%. The peptide will be identified as a weak binder if the % Rank is above the threshold of the strong binders but below the specified threshold for the weak binders, by default 2%.
   - **Novelty**: Binary value for indicating if the epitope is novel (1) or exists in the reference proteome (0). _Only present if -m flag is set to perform peptide matching in postprocessing_.
   

| Sample |  R1 |  R2 |  R3 |  Line |  chr |  allelepos |  ref |  alt |  GeneName:RefSeqID | Expression |  pos |  hla |  peptide |  core |  Of |  Gp |  Gl |  Ip |  Il |  Icore |  Identity |  Score |  Rank |  Candidate | BindLevel | Novelty |
| --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- | --- |
| test1 | 0 | 1 | 0 | line16 | chr1 | 153914523 | G | C | DENND4B:NM_014856 | 214 | 3 | HLA-B*40:01 | SERQAGAL | SERQAG-AL | 0 | 0 | 0 | 6 | 1 | SERQAGAL | line16_NM_01485 | 0.33670  | 1.30 | <= | WB |  1 |
| test1 | 1 | 1 | 0 | line8 | chr1 | 53608000 | C | T | SLC1A7:NM_001287597,SLC1A7:NM_001287595,SLC1A7:NM_006671,SLC1A7:NM_001287596 | 1560 | 2 | HLA-C*06:02 | LGFFLRTRHL | LFFLRTRHL | 0 | 1 | 1 | 0 | 0 | LGFFLRTRHL | line8_NM_001287 | 0.24655 | 1.20 | <= | WB | 1 |
| test2 | 1 | 0 | 0 | line34 | chr1 | 248402593 | C | A | OR2M4:NM_017504 | 0 | 6 | HLA-C*01:02 | VMAYERYVAI | VAYERYVAI | 0 | 1 | 1 | 0 | 0 | VMAYERYVAI | line34_NM_01750 | 0.14917 | 1.50 | <= | WB | 1 |
| test2 | 1 | 1 | 0 | line51 | chr2 | 240982213 | C | G | PRR21:NM_001080835 | NA | 2 | HLA-C*01:02 | FTHGPSSTPL | FTHPSSTPL | 0 | 3 | 1 | 0 | 0 | FTHGPSSTPL | line51_NM_00108 | 0.22570 | 0.40 | <= | SB | 1 |
| test2 | 1 | 1 | 0 | line51 | chr2 | 240982213 | C | G | PRR21:NM_001080835 | NA | 7 | HLA-C*01:02 | SSTPLHPCPF | STPLHPCPF | 0 | 1 | 1 | 0 | 0 | SSTPLHPCPF | line51_NM_00108 | 0.13137 | 2.00 | <= | WB | 1 |

   Note that if MHC-II prediction is performed, the epitope-specific columns (after **GeneName:RefID** or **Expression**) will correspond to that output, as specified [here](http://www.cbs.dtu.dk/services/NetMHCIIpan/output.php).


2. If there are not multiple regions from a single patient the resulting summary table will appear as follows (the following are the same for both multiregion below and single region):
   - **Sample**: Sample identifier
   - **Total**: Total Neoantigen burdens that are of proper range.
   - **Total_WB**: Total Neoantigen burdens of weak binding affinity.
   - **Total_SB**: Total Neoantigen burdens of strong binding affinity.

| Sample | Total | Total_WB | Total_SB |
|  --- |  --- |  --- |  --- |
|  Pat1 |  72 |  72 |  0 |
|  Pat2 |  33 |  23 |  10 |

3. If multiple regions are specified then the output will look as follows (scroll left or right to view all):
   - For cases of multiregion samples, the same information for totals are given, but also for each region in the vcf.
   - Heterogeneity (e.g. clonal, subclonal, and shared) information is also measured and printed out. This yields counts of clonal subclonal and shared.
     - For _shared_ neoantigens there must be >2 regions present, otherwise shared will be 0. **This pipeline can handle samples with different numbers of regions**.

| Sample | Total | Total_WB | Total_SB | Total_Region_1 | Total_Region_n | Total_WB_Region_1 | Total_WB_Region_n | Total_SB_Region_1 | Total_SB_Region_n | Clonal | Subclonal | Shared | Clonal_WB | Clonal_SB | Subclonal_WB | Subclonal_SB | Shared_WB | Shared_SB |
|  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |  --- |
| test1 | 86 | 65 | 21 | 48 | 51 | 0 | 36 | 40 | 0 | 12 | 11 | 0 | 13 | 73 | 0 | 11 | 2 | 54 | 19 | 0 | 0 |
| test2 | 86 | 66 | 20 | 57 | 43 | 0 | 46 | 30 | 0 | 11 | 13 | 0 | 14 | 72 | 0 | 10 | 4 | 56 | 16 | 0 | 0 |

4. The above two files are reported separately for single nucleotide changes and indels (and/or other genetic alterations resulting in more than 1 amino acid change).
   - _ExperimentName_.neoantigens.txt and _ExperimentName_.neoantigens.summarytable.txt contain single amino acid changes.
   - _ExperimentName_.neoantigens.Indels.txt and _ExperimentName_.neoantigens.Indels.summarytable.txt contain neoantigen information arising from indel/frameshift/stop-loss events.

