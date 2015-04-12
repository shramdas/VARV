# VARV

Authors:

* Shweta Ramdas (sramdas@umich.edu)
* Ryan Welch (welchr@umich.edu)
* Tanya Teslovich (teslo@umich.edu) 

Questions should be sent to all three. Bugs/problems can be posted to the issues section of this repository. 

## Contents

* [Synopsis](#synopsis)
* [Installation](#installation)
* [Usage](#usage)
* [Basic Configuration File](#configuration)
* [Additional Options](#additional-options)
* [Example](#example)
* [File formats](#file-formats)
* [Output](#output)
* [Plots](#plots)
* [Options](#options)
* [Limitations](#limitations)
* [License](#license)

## Synopsis

VARV is a pipeline for automating many of the steps required for performing gene or group based tests with [EPACTS](http://genome.sph.umich.edu/wiki/EPACTS). VARV is a command line tool that can be run from a Linux system.

## Installation

VARV requires: 

* [EPACTS](http://genome.sph.umich.edu/wiki/EPACTS)
* Python 2.7.x (untested on / likely does not work with 3.x)
* R 3.1 or greater
* Tabix (available as a part of SAMtools, http://samtools.sourceforge.net/) located somewhere on your $PATH
* An internet connection during setup
* An up to date version of Chrome or Firefox (for viewing the output report/visualizations)

To download, grab the tarball from the [Releases](https://github.com/shramdas/varv/releases) page. 

To install, simply extract the tarball, and then execute the setup script: 

```bash
tar zxf varv.tar.gz 
cd varv/bin/
python setup.py 
```

The setup script will attempt to install a number of required python packages into a virtualenv located in `varv/env`. This requires internet access on port 443 for SSL traffic. 

## Usage 

The main program is located in `bin/varv` under the extracted directory. To run this without needing to directly specify it, you can create a symbolic link to it from somewhere accessible on your `PATH`, usually something like: 

```
sudo ln -s -t /usr/local/bin /path/to/varv_dir/bin/varv
```

Then you can invoke the program on a configuration file like so: 

```
varv your_config.cfg
```

## Configuration File

The config file describes the input parameters and files to be submitted to VARV. The config file is in the METAL format, with each line containing a parameter and value. 
The necessary parameters of the config file are:

```
INPUTDIR  The input directory containing all the input files.

VCFFILE The input vcf file (thie file does not have to be in the INPUTDIR). Only the file name has to be specified if file is in INPUTDIR (see additional options), otherwise the absolute path. 

PEDFILE File containing phenotype information for samples (tab-delimited). Assumed to be in INPUTDIR.

MODEL The model for analysis in the format 'phenotype ~ covariate2 + covariate2' Ex: phenotype ~ sex, or only phenotype (if no covariate is to be used.)

TEST  Which test is to be used for group-wise test. Any of the options for EPACTS can be specified. If multiple tests are to be run, use multiple TEST lines (Ex: TEST  group=skat)
OUTPREFIX Output prefix for this run.

PROCESS Execute the set of instruction specified above this line (and after the previous instance of PROCESS)

```

## Additional Options

```

EPACTSDIR Location of EPACTS if EPACTS is not on system path

KINSHIPFILE Kinship file, if one has been generated previously

ANNOTFILE Annotation for variants, one variant per line. 

PVALUETHRESHOLD 

GENELIST

GROUPFILE

SINGLEMARKERTEST

```

## Example

Below is an example of a config file, along with comments describing each of the options. You can remove the comments for brevity in your actual config file. 

The config file is whitespace delimited, so the arguments to each option must not contain spaces. 

```
# Location of tabix. Can just be "tabix" if it exists on your path. 
TABIX tabix

# Location of EPACTS. 
EPACTSDIR	/net/fantasia/home/hmkang/bin/epactsRelease/bin

# Directory containing input files (VCFs, or any other file 
# you wish to refer to simply by name and not full path.) 
INPUTDIR	/net/snowwhite/home/xlsim/METSIM/FFA/gene-level/gene_plots

# Either absolute path to a VCF file, or just the name if it 
# is located in the INPUTDIR. 
VCFFILE	/net/snowwhite/home/teslo/ExomeChip/Gene-level_analysis/Tanya_dataset/omni_exchip_GoT2D_mean-genotype_10042_FINAL.vcf.gz

# PED file containing phenotype information. 
PEDFILE	/net/snowwhite/home/welchr/projects/FFA/pheno/excluding_all_t2d/ffa_phenos_studyage-1_BMI-1_smokefix-1_dropstatin-0_adjuststatin-1_dropt2d-1.ped

# Kinship file, if one has already been generated previously. Otherwise it will be created automatically 
# if you select a test that requires a kinship file. 
KINSHIPFILE	/net/snowwhite/home/aujackso/METSIM/OMNI_Express/Geno/metsim_omni.kinf

# Annotations for variants. This can be used for automatically generating a group/mask file, or for showing annotations 
# in resulting output files. 
ANNOTFILE /net/snowwhite/home/welchr/projects/varv/tests/omni_exchip.hg19.vep.annotation.gz

# Which columns from the annotation file do we care about showing in the output? 
# This can be omitted, in which case, no annotations will be shown for variants in the output. 
ANNOTCOLUMNS	vepFUNC

# Specify a model here for analysis. The usual R type of formula should be used. 
# For example: PC ~ AGE + BMI, or just PC if no covariates are needed. 
MODEL	PC

# Gene-based test threshold for significance. 
# Set to 1 to see everything, but this may cause the plotting output to be overwhelming. 
# If not specified, the default is 1e-4. 
PVALUETHRESHOLD	1

# If you wish to restrict analysis to only a few genes, you can list them here. 
# Note that you'll likely want to set PVALUETHRESHOLD to 1 (as above) to be sure
# you get output for all of the genes listed, regardless of their significance. 
GENELIST LIPG,LCAT,FADS1

# Output directory + prefix for output files. 
OUTPREFIX	/net/snowwhite/home/welchr/projects/varv/tests/testcfg/blah

# Group file specifying the mapping from group names (genes) to sets of variants. 
# If this is omitted, it will be generated automatically. 
GROUPFILE PTV+missense.01.grp

# The pipeline will also execute single marker tests, to include them with the results of the
# gene-based tests for interpretation. This specifies the name of the test as EPACTS would expect it. 
# http://genome.sph.umich.edu/wiki/EPACTS#Single_Variant_Tests
SINGLEMARKERTEST	q.emmax

# Specify the group based test that should be run. If you want to execute multiple tests, 
# just use multiple TEST lines. 
# http://genome.sph.umich.edu/wiki/EPACTS#Gene-wise_or_group-wise_tests
TEST group=mmskat
TEST group=VT

# If you specify a SKAT-based test, and you want to use the SKAT-O optimal test, 
# use this option. 
SKATO	ON

# Parallel options for our local MOSIX cluster. This is likely not useful to anyone else!
#NJOBS 100
#MOSIX --mosix-nodes="\`/net/fantasia/home/gjun/bin/pick_mainnode\`"

# Execute this set of instructions. 
PROCESS

# Another run. 
# All of the parameters set are remembered, so if you wish to run another analysis with only one 
# or few parameters changed, just change them and issue PROCESS again (see below for an example.) 
# The exception is that TEST must be re-specified to be explicit about which tests should be run. 
MODEL FFA ~ BMI + AGE + PC1 + PC2
OUTPREFIX /my/new/analysis/prefix
TEST group=mmskat

# If you want to reset an option back to its default: 
PVALUETHRESHOLD RESET

PROCESS

```

The following additional config file options exist: 

```

# If you've provided an annotation file, but want to perform a filter on it, use
# this option. The filter is written with a simple python syntax. For more details, 
# see: http://pandas.pydata.org/pandas-docs/version/0.15.1/indexing.html#indexing-query
FILTERANNOT (vepGENE == 'LIPG') & (CAD > 0.4)

# Is there an alternative genotype field that EPACTS should use? 
FIELD EC

# If VCF file is split up into chromosomes, specify this option. 
SEPCHR ON

# Only consider genes with a total MAC > this value. 
GENEMINMAC 5

# Only consider genes with > this number of variants. 
MINVARS 2

# Forcefully insert a string of command line arguments into the EPACTS command line
# when running group/gene tests
EPACTSCMD --some-opt 0.1

# Place restrictions on min/max MAF and MAC to be used during both
# single variant and group based association tests. 
# These are the current defaults: 
MIN_MAC 1
MIN_MAF 0
MAX_MAF 1

# Do you want to see more verbose output? This will spew every EPACTS command line out, 
# along with some other (potentially) useful debug output. 
VERBOSE ON

# REML file as previously generated by EPACTS. 
# This should only be used with caution, and only when this *exact same* 
# analysis has previously been executed. It is dependent on just about 
# everything - genotypes, phenotype, covariates, kinship. 
REMLFILE /net/snowwhite/home/xlsim/METSIM/FFA/gene-level/results_noT2D_adjmed/emmaxSKATO/PC_PTV+missense.01.grp_emmaxSKATO.reml

```

## File formats

### PEDFILE

[PED file for phenotypes and covariates](http://genome.sph.umich.edu/wiki/EPACTS#PED_file_for_Phenotypes_and_Covariates)

### VCFFILE

[VCF format 4.1 or higher](http://samtools.github.io/hts-specs/)

### GROUPFILE

[Creating a marker group file](http://genome.sph.umich.edu/wiki/EPACTS#Creating_marker_group_file)

### ANNOTFILE

Basically just the first few columns of a VCF file, along with additional columns specifying the annotations: 

| CHROM | POS    | ID         | REF | ALT | QUAL | FILTER | INFO                                   | FORMAT | ALLELEA | ALLELEB | TOP | vepFUNC                 | vepGENE        | vepNS | vepPTV | CADD  |
|-------|--------|------------|-----|-----|------|--------|----------------------------------------|--------|---------|---------|-----|-------------------------|----------------|-------|--------|-------|
| 1     | 752566 | rs3094315  | G   | A   | 100  | PASS   | NS=10042;AC=16000;AN=20084;AF=0.796654 | GT     | A       | G       | +   | intron_variant          | RP11-206L10.10 | 0     | 0      | 0.189 |
| 1     | 752721 | rs3131972  | A   | G   | 100  | PASS   | NS=10042;AC=15994;AN=20084;AF=0.796355 | GT     | A       | G       | +   | intron_variant          | RP11-206L10.10 | 0     | 0      | 3.952 |
| 1     | 762320 | exm2268640 | C   | T   | 100  | PASS   | NS=10042;AC=6;AN=20084;AF=0.000299     | GT     | T       | C       | -   | non_coding_exon_variant | LINC00115      | 0     | 0      | 3.28  |
| 1     | 776546 | rs12124819 | A   | G   | 100  | PASS   | NS=10042;AC=4295;AN=20084;AF=0.213852  | GT     | A       | G       | +   | intron_variant          | RP11-206L10.11 | 0     | 0      | 1.62  |
| 1     | 798959 | rs11240777 | G   | A   | 100  | PASS   | NS=10042;AC=5484;AN=20084;AF=0.273053  | GT     | A       | G       | +   | downstream_gene_variant | RP11-206L10.11 | 0     | 0      | 0.265 |
| 1     | 800007 | rs6681049  | T   | C   | 100  | PASS   | NS=10042;AC=20052;AN=20084;AF=0.998407 | GT     | T       | C       | -   | downstream_gene_variant | FAM41C         | 0     | 0      | 0.601 |
| 1     | 838555 | rs4970383  | C   | A   | 100  | PASS   | NS=10042;AC=6335;AN=20084;AF=0.315425  | GT     | A       | C       | +   | upstream_gene_variant   | RP11-54O7.16   | 0     | 0      | 1.012 |
| 1     | 846808 | rs4475691  | C   | T   | 100  | PASS   | NS=10042;AC=3658;AN=20084;AF=0.182135  | GT     | T       | C       | -   | intron_variant          | RP11-54O7.16   | 0     | 0      | 5.381 |
| 1     | 854250 | rs7537756  | A   | G   | 100  | PASS   | NS=10042;AC=3657;AN=20084;AF=0.182085  | GT     | A       | G       | +   | non_coding_exon_variant | RP11-54O7.3    | 0     | 0      | 1.108 |

### KINSHIPFILE

This is automatically created by the pipeline if a test requiring a kinship matrix is requested. However, if you wish to make your own kinship file, see [Creating the kinship matrix](http://genome.sph.umich.edu/wiki/EPACTS).

## Output

The pipeline creates a number of output files that are used during the analysis, and may be potentially useful for interpreting the results. 

`prefix.singlemarker.epacts.gz`

The single marker association results (for those variants going into gene based tests)

`prefix.groupfile.filtered.txt`

The final group file entering into the gene based testing after filtering has (potentially) been applied. For example, variants may be excluded for failing to meet `MAX_MAF`. 

`prefix.genotype_trait_table.gene-?.tab` 

For each significant gene (according to `PVALUETHRESHOLD`), a file will be created containing the genotypes for the variants going into that gene's test, along with the trait values for each individual, like so: 

| IND_ID   | PC           | 18:47091686_G/A | 18:47095862_C/T | 18:47101838_G/A |
|----------|--------------|-----------------|-----------------|-----------------|
| 10741    | -3.85229944  | 0/0             | 0/0             | 0/0             |
| MJQQ-901 | -3.574392931 | 0/0             | 0/0             | 0/0             |
| MCMP-901 | -3.438443305 | 0/0             | 0/0             | 0/0             |
| MFID-901 | -3.346259618 | 0/0             | 0/0             | 0/0             |
| 8071     | -3.275935989 | 0/0             | 0/0             | 0/0             |
| MJNO-901 | -3.218829067 | 0/0             | 0/0             | 0/0             |
| MJEG-901 | -3.170612436 | 0/0             | 0/0             | 0/0             |
| 11618    | -3.128801711 | 0/0             | 0/0             | 0/0             |
| 11577    | -3.09183487  | 0/0             | 0/0             | 0/0             |

`prefix.pheno.ped`

The phenotype file that entered into both the single variant and gene based tests, after any potential filtering was applied.

`prefix.plots/` 

Directory containing files for creating some of the QC plots and visualizations. In particular: 

* `top_genes.html` contains a report showing QQ plots, manhattan plots, and a graphic to help interpret the gene-based tests that shows each individual carrier of a rare allele for a variant contributing to a significant gene test (see [Plots](#plots) for more info.) 

* `plot_variants.tsv` has the genotype data melted into long format for use in the `top_genes.html` figure and looks like the following: 

| VARIANT         | IND     | GENOTYPE | RARE_COUNT | TRAIT  |
|-----------------|---------|----------|------------|--------|
| 16:67976823_C/T | SAMPLE1 | 0/1      | 1          | 1.69   |
| 16:67976823_C/T | SAMPLE2 | 0/1      | 1          | -0.138 |
| 16:67976823_C/T | SAMPLE3 | 0/1      | 1          | -0.409 |
| 16:67976823_C/T | SAMPLE4 | 0/1      | 1          | 0.16   |

* `plot_genes.tsv` contains the table exactly as seen in `top_genes.html`, but in tab-delimited format. For example: 

| GENE | mmskat.P | VARIANT         | SV.P     | BETA    | MAF      | MAC | vepFUNC                 |
|------|----------|-----------------|----------|---------|----------|-----|-------------------------|
| LIPG | 6.20E-11 | 18:47091686_G/A | 2.96E-06 | 0.76    | 0.00225  | 38  | splice_acceptor_variant |
| LIPG | 6.20E-11 | 18:47095862_C/T | 0.000184 | 0.437   | 0.0045   | 76  | missense_variant        |
| LIPG | 6.20E-11 | 18:47101838_G/A | 0.792    | -0.263  | 6.00E-05 | 1   | missense_variant        |
| LIPG | 6.20E-11 | 18:47109939_G/A | 0.957    | -0.0192 | 0.00047  | 8   | missense_variant        |
| LIPG | 6.20E-11 | 18:47109955_A/G | 4.36E-05 | 0.37    | 0.00728  | 123 | missense_variant        |

`prefix.?.epacts`

An EPACTS file for the gene-based test results, one per test requested. 

## Plots

The `top_genes.html` file contains 3 different visualizations for QC and interpretation of your results. 

* Carriers of rare alleles for single variants contributing to significant gene-based tests

![carrier-plot](http://www.umich.edu/~welchr/top_genes_gene_and_single.png)

Each row represents a variant contributing to a significant gene-based test (SKAT P=6.2e-11 for LIPG.) The symbols in the leftmost column of the table correspond to each individual carrier of that variant, and are aligned with the x-axis of the histogram above to show where each variant carrier lies in the phenotype distribution. An interactive tooltip displayed over each carrier shows the phenotype value and the sample identifier. 

* QQ and Manhattan plots

These plots are created for two subsets of variants: 

1. Variants assigned to genes that pass our inclusion filters (such as minimum total MAC or minimum number of variants)
2. Variants assigned to genes failing filters

The manhattan plots show single variant association statistics. 

![qqmanhat-plots](http://www.umich.edu/~welchr/top_genes_qq_manhat.png)

## License

Copyright (C) 2014 Shweta Ramdas, Ryan Welch, The University of Michigan

VARV is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

VARV is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
