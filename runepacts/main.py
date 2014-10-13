#!/usr/bin/env python

#===============================================================================
# Copyright (C) 2013 Shweta Ramdas, Ryan Welch, The University of Michigan
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================


def bootstrap_lib(lib_path):
	import os, sys
	from glob import glob

	exist_eggs = [os.path.basename(i) for i in sys.path];

	sys.path.insert(1,lib_path);

	for d in glob(os.path.join(lib_path,"*egg")):
		if os.path.basename(d) in exist_eggs:
			continue;
		else:
			sys.path.insert(1,d);

bootstrap_lib("/net/snowwhite/home/welchr/lib/python2.7/site-packages");

import os, sys, re, getopt, subprocess, gzip, numpy, time, optparse, logging, pprint, signal
import pandas
import pandas.computation.ops
import runepacts
from copy import deepcopy
from runepacts.util import *
from distutils.dir_util import mkpath
import functools as ft

_RUNEPACTS_DEBUG = True

EPACTS_KIN_FILTER_ARGS = "-min-maf 0.01 -min-callrate 0.95"
EPACTS_MAKE_GROUP_ARGS = "-nonsyn"
KINSHIP_TESTS = "q.emmax emmax emmaxVT emmaxCMC emmaxSKAT mmskat".split()

# Freeze the verbose argument of run_bash to be dependent on whether
# we're in debug mode.
# Now you can just call run_bash(cmd) without specifying verbose.
run_bash = ft.partial(run_bash,verbose=_RUNEPACTS_DEBUG)

def get_defaults():
	"""
	Get a dictionary of default values for program options.
	"""

	defaults = dict()

	defaults['PEDCOLUMNS'] = []
	defaults['PVALUETHRESHOLD'] = 0.05
	defaults['SEPCHR'] = False
	defaults['SINGLEMARKERTEST'] = 'q.linear'
	defaults['FILTERMAF'] = 0.05
	defaults['MINMAF'] = 0
	defaults['GENEMINMAC'] = 5
	defaults['MAXMAF'] = 1
	defaults['VERBOSE'] = 'OFF'
	defaults['GENELIST'] = None
	defaults['MINMAC'] = 0
	# defaults['ANNOTGENECOL'] = 4
	# defaults['ANNOTVARCOL'] = 0
	# defaults['ANNOTPOSCOL'] = 1
	defaults['ANNOT_GENE_COL'] = "vepGENE"
	defaults['ANNOT_CHROM_COL'] = "CHROM"
	defaults['ANNOT_POS_COL'] = "POS"
	defaults['ANNOT_REF_COL'] = "REF"
	defaults['ANNOT_ALT_COL'] = "ALT"
	defaults['MARKERMINMAF'] = 0
	defaults['MARKERMAXMAF'] = 1
	defaults['MARKERMINMAC'] = 1
	defaults['MINVARS'] = 2
	defaults['EPACTS'] = 'epacts'
	defaults['EPACTSCMD'] = ''
	defaults['NJOBS'] = 1
	defaults['SKATO'] = False

	#defaults['MOSIX'] = '--mosix-nodes="\`/net/fantasia/home/gjun/bin/pick_mainnode\`"'

	return defaults

def run_plots(options,phenotype,covariates,tests,pvalt):
	log_key = gen_analysis_key(options)
	info = logging.getLogger(log_key).info
	warn = logging.getLogger(log_key).warning
	debug = logging.getLogger(log_key).debug

	# These files must have been previously created by the pipeline.
	required_for_plot = map(lambda x: options["OUTPREFIX"] + x,[".allvariants.txt",".pheno.ped",".summary.txt",".genes_passing_filters.txt"]);

	if all(map(os.path.isfile,required_for_plot)):
		plot_loc = os.path.join(runepacts.__path__[0],"R/plots.R");
		if len(covariates) > 0:
			cmd = "R --vanilla --slave --args NULL {phenotype} {out_prefix} {tests} {pval_thresh} {covars} < {script}";
			run_cmd = cmd.format(
				out_prefix = options["OUTPREFIX"],
				phenotype = phenotype,
				tests = tests,
				pval_thresh = str(pvalt),
				covars = options['MODEL'].split('~')[1],
				script = plot_loc
			)

			run_bash(run_cmd)

		else:
			cmd = "R --vanilla --slave --args NULL {phenotype} {out_prefix} {tests} {pval_thresh} NA < {script}"
			run_cmd = cmd.format(
				out_prefix = options["OUTPREFIX"],
				phenotype = phenotype,
				tests = tests,
				pval_thresh = str(pvalt),
				script = plot_loc
			)

			run_bash(run_cmd)

		if _RUNEPACTS_DEBUG:
			debug(run_cmd)

		info("Output written to: " + options['OUTPREFIX'] + '.topgenes.pdf')

	else:
		warn("Skipping plot, could not find required files: \n%s" % "\n".join(required_for_plot))

def check_analysis(adict):
	"""
	Check an analysis for validity.
	Raises an exception if it is not valid.

	:param adict: Dictionary with options for running this analysis.
	:raises: ValueError
	"""

	# Required options
	for k in "MODEL OUTPREFIX TEST".split():
		if k not in adict:
			raise ValueError, "Error: option '%s' must be set in your config file" % k

	if len(adict["TEST"]) == 0:
		raise ValueError, "Error: must specify at least 1 TEST for each PROCESS block in your config file."

def read_config(filepath):
	"""
	Reads a config file into a list of dictionaries, one dictionary per analyses containing setting/value pairs.
	:param filepath: Full path to the config file
	:return: A list of dictionaries
	"""

	defaults = get_defaults()
	analyses = []

	current = deepcopy(defaults)
	with open(filepath) as f:
		for line in f:
			# Skip blank lines
			if line.strip() == "":
				continue

			# Skip commented lines.
			if line.startswith("#"):
				continue

			# We now have a full analysis, save it
			if line.startswith("PROCESS"):
				# Check it for validity first.
				check_analysis(current)

				# It's valid, so keep it.
				analyses.append(current)

				# Copy it, since all previous settings are still valid, and further lines can modify them.
				current = deepcopy(current)

				# Reset the tests to be done. This is apparently supposed to be specified within each PROCESS block.
				current["TEST"] = []

				continue

			# Check length (should get two elements, one key and one value)
			ls = line.strip().split()
			if len(ls) < 2:
				print >> sys.stderr, "Ignoring line in config, expected two elements, got only one: %s" % line.rstrip()
				continue

			# This splits on any whitespace, but only does 1 split. The first leftmost split is done.
			# So for example:
			# MODEL      phenotype ~ covar1 + covar2
			# Turns into:
			# ["MODEL","phenotype ~ covar1 + covar2"]
			key, value = line.split(None,1)
			value = value.strip()

			# User can specify TEST multiple times, one per line, denoting each test to be performed
			if key == "TEST":
				current.setdefault(key,[]).append(value)
				continue

			if key == "MODEL":
				if "~" in value:
					covars = value.split("~")[1].split("+")
					covars = map(str.strip,covars)
					pheno = value.split("~")[0].strip()
				else:
					covars = []
					pheno = value.strip()

				current["COVARIATES"] = covars
				current["PHENOTYPE"] = pheno
				current["MODEL"] = value

				continue

			if key == "INPUTDIR":
				if not os.path.isdir(value):
					die("Error: input directory does not exist: %s" % value)

				# Assume this is the VCFDIR as well, unless we later encounter VCFFILE
				# and it contains directories in its path
				current["INPUTDIR"] = value
				current["VCFDIR"] = value

				continue

			if key == "VCFFILE":
				if '/' in value or '\\' in value:
					check_file_exists(value)

					vcf_dir, vcf_file = os.path.split(value)
					current["VCFDIR"] = vcf_dir
					current["VCFFILE"] = value
				else:
					current["VCFDIR"] = current["INPUTDIR"]
					current["VCFFILE"] = os.path.join(current["INPUTDIR"],value)

					check_file_exists(current["VCFFILE"])

				continue

			if key in "GROUPFILE REMLFILE KINSHIPFILE ANNOTFILE PEDFILE".split():
				if '/' in value or '\\' in value:
					current[key] = value
				else:
					current[key] = os.path.join(current["INPUTDIR"],value)

				check_file_exists(current[key])

				continue

			if key == "GENELIST":
				if os.path.isfile(value):
					current[key] = value
				elif "/" in value or "\\" in value:
					current[key] = os.path.join(current["INPUTDIR"],value)
				else:
					current[key] = value

				continue

			# This is apparently a "reset" back to default
			if value in ("NA","NaN","None","."):
				value = defaults[key]

			# Handle specific cases
			if key in "ANNOTGENECOL ANNOTVARCOL ANNOTPOSCOL".split():
				value = int(value) - 1

			if key == "MARKERMINMAF":
				value = float(value)

			if key == "MARKERMINMAC":
				value = int(value)
				if value < 0:
					die("Error: invalid MARKERMINMAC in config file: %s" % str(value))

			if key == "GENEMINMAC":
				value = int(value)
				if value < 0:
					die("Error: invalid GENEMINMAC: %s" % str(value))

			if key == "MINVARS":
				value = int(value)
				if value < 0:
					die("Error: invalid MINVARS: %s" % str(value))

			if key == "EPACTSCMD":
				value = " %s " % value.strip("'")

			if key == "EPACTSDIR":
				current["EPACTS"] = os.path.join(value,"epacts")

			if key == "SEPCHR":
				value = True

			if key == "SKATO":
				value = True

			if key == "PVALUETHRESHOLD":
				value = float(value)
				if value > 1 or value < 0:
					die("Error: invalid p-value threshold specified in config file: " % str(value))

			if key == "PEDCOLUMNS":
				value = value.split(",")

			current[key] = value

	return analyses

def gen_analysis_key(an):
	import hashlib

	return hashlib.md5("".join([str(x) for x in an.itervalues()])).hexdigest()

def main():
	op = optparse.OptionParser()
	op.add_option("--plotonly",help="Only create the plot PDF, assuming you've already run the pipeline once.",default=False,action="store_true")

	opts, args = op.parse_args()

	if len(args) == 0:
		die("Error: must specify config file, e.g. runepacts.py config.cfg")

	if not os.path.isfile(args[0]):
		die("Error: config file does not exist: %s" % args[0])

	# Read in the config file into a list of dictionaries, each one specifying a particular analysis to run
	analyses = read_config(args[0])

	# Loop over each analysis specified in the config file. An analysis is just a set of instructions between PROCESS
	# commands.
	for aopts in analyses:
		tests = aopts["TEST"]

		if len(tests) == 0:
			die("Error: must specify at least 1 test for each PROCESS.")

		outprefix_dir = os.path.split(aopts["OUTPREFIX"])[0]
		if not os.path.isdir(outprefix_dir):
			mkpath(outprefix_dir)

		# A unique key for this particular analysis. Just makes it easier to setup a logger.
		analysis_key = gen_analysis_key(aopts)

		# Create a logger to log to both STDOUT and a log file
		logger = logging.getLogger(analysis_key)

		# Setup logging to file
		ffmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
		fhandler = logging.FileHandler(aopts["OUTPREFIX"] + ".log")
		fhandler.setFormatter(ffmt)
		logger.addHandler(fhandler)

		# Setup logging to STDOUT
		logger.addHandler(logging.StreamHandler(sys.stdout))

		# If we're in debug mode, print everything.
		if _RUNEPACTS_DEBUG:
			logger.setLevel(logging.DEBUG)
		else:
			logger.setLevel(logging.INFO)

		logger.info("Running analysis: ")
		logger.info(pprint.pformat(aopts))

		orig_vcf_path = aopts["VCFFILE"]
		epacts = aopts["EPACTS"]
		pedcolumns = aopts["PEDCOLUMNS"]
		user_cmd = aopts.get("EPACTSCMD","")
		field = '-field %s' % aopts["FIELD"] if "FIELD" in aopts else ' '
		min_maf = " -min-maf %s" % str(aopts["MINMAF"])
		marker_min_maf = '-min-maf ' + str(aopts['MARKERMINMAF'])
		marker_min_mac = '-min-mac ' + str(aopts['MARKERMINMAC'])
		marker_max_maf = '-max-maf ' + str(aopts['MARKERMAXMAF'])
		covariates = aopts["COVARIATES"]
		phenotype = aopts["PHENOTYPE"]

		tests_to_write = ",".join([t.split("=")[1] for t in tests])

		if opts.plotonly:
			logger.info("Creating plots...")
			run_plots(aopts,phenotype,covariates,tests_to_write,aopts["PVALUETHRESHOLD"])

			continue

		# First step: Index VCF if tabix index file does not already exist
		# if sepchr is on, then tabix all other chromosomes too
		if aopts["SEPCHR"]:
			sep_vcfs = match_sepchr_vcfs(orig_vcf_path)
			for f in sep_vcfs:
				if os.path.isfile(f + ".tbi"):
					continue

				tabix_cmd = "tabix -p vcf %s" % f
				logger.debug(tabix_cmd)
				run_bash(tabix_cmd)

		elif not os.path.isfile(orig_vcf_path + '.tbi'):
			logger.info("Indexing VCF...")

			tabix_command = 'tabix -p vcf -f ' + orig_vcf_path
			logger.debug(tabix_command)
			run_bash(tabix_command)

		# The VCF used for tests may end up being different than the original full VCF file
		# For example, if filtered down to a specific set of genes using GENELIST
		vcf_for_tests = orig_vcf_path

		# Then, make a ped file after using the filters provided
		logger.info("Creating PED file...")

		pedfile = pandas.read_table(aopts['PEDFILE'],sep="\t",header=0)

		# The variable filter contains the awk filter to use for the ped file. we are replacing the filter-column name with the column
		# number to get the awk command
		logger.info("Filtering PED file...")

		if 'FILTERPED' in aopts:
			pfilter = aopts['FILTERPED']

			try:
				pedfile = pedfile.query(pfilter)
			except pandas.computation.ops.UndefinedVariableError:
				die("Error: PED file filter specified a variable that does not exist: %s" % pfilter,logger)

			if pedfile is None or pedfile.shape[0] == 0:
				# Filter resulted in an empty PED file, we can't continue.
				die("Filter '%s' on PED file resulted in an empty data set" % pfilter,logger)

		pedfile.iloc[:,1].to_csv(aopts['OUTPREFIX'] + ".samples_to_keep.txt", sep="\t",index=False, na_rep='NA')
		keep_samples = pedfile.iloc[:,1]
		if len(pedcolumns) > 0:
			if len(covariates) > 0:
				for i in range(0,len(covariates)):
					if covariates[i].strip() not in pedcolumns:
						pedcolumns.append(covariates[i].strip())

			pedfile = pedfile.loc[:,pedcolumns]

		# Make sure PED columns are named correctly.
		pedfile.columns = ["#FAM_ID","IND_ID","FAT_ID","MOT_ID","SEX"] + list(pedfile.columns[5:])

		# Checking if there are categorical covariates. If yes, convert to dummy variables
		logger.info("Adding covariates...")
		covar_cmd = ''

		for i in range(0,len(covariates)):
			# if covariate is categorical, then create dummy variables for it
			if 'as.factor' in covariates[i]:
				covariate = covariates[i].replace('as.factor(', '')
				covariate.strip(')')
				values = list(set(list(pedfile[covariate])))			# get the categories
				numvalues = len(set(list(pedfile[covariate])))			# get the number of categories
				covcolumns = numpy.empty([len(pedfile),numvalues-1])

				for row in range(0,len(pedfile)):
					for col in range(0,numvalues-1):
						if pedfile[covariate][row] == values[col]:
							covcolumns[row][col] = 1
						else:
							covcolumns[row][col] = 0

				num = 0
				for _ in range(0,len(values)-1):
					pedfile[covariate + '_' + str(num)] = covcolumns[:,num]
					num += 1
					covar_cmd = covar_cmd + ' -cov ' + covariate + '_' + str(num)

			else:		# no dummy variable has to be added to the file
				covar_cmd = covar_cmd + ' -cov ' + covariates[i] + ' '

		pedfile.to_csv(aopts['OUTPREFIX'] + '.pheno.ped', sep="\t",index=False, na_rep='NA')

		# Create groupfile if one is not already given by the user
		if 'GROUPFILE' not in aopts:
			logger.info("Creating groupfile...")

			final_groupfile_name = aopts['OUTPREFIX'] + '.groupfile.txt'

			# No annotation file provided by the user, so we have to use epacts to create the groupfile
			if 'ANNOTFILE' not in aopts:
				logger.info("No annotation file provided, using EPACTS to create group file")

				if not aopts["SEPCHR"]:
					annot_cmd = "{epacts} anno -in {vcf} -out {out}".format(
						epacts = epacts,
						vcf = orig_vcf_path,
						out = aopts['OUTPREFIX'] + '.anno.vcf.gz'
					)

					logger.debug(annot_cmd)
					run_bash(annot_cmd)

					annot_cmd = "{epacts} make-group -vcf {vcf} -out {out} --format epacts {args}".format(
						epacts = epacts,
						vcf = aopts['OUTPREFIX'] + '.anno.vcf.gz',
						out = final_groupfile_name,
						args = EPACTS_MAKE_GROUP_ARGS
					)
					logger.debug(annot_cmd)
					run_bash(annot_cmd)

				# SEPCHR, so we annotate each chromosome separately
				else:
					sep_vcfs = match_sepchr_vcfs(aopts["VCFFILE"])

					try:
						os.remove(final_groupfile_name)
					except:
						pass

					# Loop over each of the VCF files separated by chromosome
					for sep_vcf_fpath in sep_vcfs:
						sep_vcf_name = os.path.split(sep_vcf_fpath)[1]

						annot_cmd = "{epacts} anno -in {vcf} -out {out}".format(
							epacts = epacts,
							vcf = sep_vcf_fpath,
							out = aopts['OUTPREFIX'] + '.anno.' + sep_vcf_name
						)
						logger.debug(annot_cmd)
						run_bash(annot_cmd)

						annot_cmd = "{epacts} make-group -vcf {vcf} -out {out} -format epacts {args}".format(
							epacts = epacts,
							vcf = aopts['OUTPREFIX'] + '.anno.' + sep_vcf_name,
							out = aopts['OUTPREFIX'] + '.group.' + sep_vcf_name,
							args = EPACTS_MAKE_GROUP_ARGS
						)
						logger.debug(annot_cmd)
						run_bash(annot_cmd)

						annot_cmd = "cat {0} >> {1}".format(
							aopts['OUTPREFIX'] + '.group.' + sep_vcf_name,
							final_groupfile_name
						)
						logger.debug(annot_cmd)
						run_bash(annot_cmd)

			# They specified an annotation file, so we should use that instead of EPACTS' anno.
			else:
				# A list of extra user-specified columns to keep along with the annotations.
				annot_extra_cols = aopts["ANNOTCOLUMNS"].split(",") if "ANNOTCOLUMNS" in aopts else None

				# Load in the annotation file.
				# If FILTERANNOT was specified, it will be used to filter the data frame.
				annofile = load_annotation(
					aopts["ANNOTFILE"],
					aopts['ANNOT_CHROM_COL'],
					aopts['ANNOT_POS_COL'],
					aopts['ANNOT_REF_COL'],
					aopts['ANNOT_ALT_COL'],
					aopts['ANNOT_GENE_COL'],
					annot_extra_cols,
					aopts.get("FILTERANNOT",None)
				)

				# Create the group file from the (possibly filtered) annotations.
				create_group_file(
					final_groupfile_name,
					annofile,
					aopts['ANNOT_GENE_COL'],
					logger = logger
				)

		else:
			logger.info("Using existing groupfile...")
			final_groupfile_name = aopts['GROUPFILE']

		# This code subsets the input VCF files down to only those genes in GENELIST, or if GENELIST wasn't specified,
		# then down to only variants within genes in the group file (either the generated one, or the one provided by the user.)
		if aopts["GENELIST"] is None:
			genes_to_keep = groupfile_get_genes(final_groupfile_name)
		else:
			genes_to_keep = []

			# If GENELIST is a file, read it in.
			# Otherwise it's a comma-delimited string of genes.
			gene_list_file = aopts['GENELIST']
			if os.path.isfile(gene_list_file):
				with open(gene_list_file) as f:
					for line in f:
						genes_to_keep.append(line.strip())
			else:
				genes_to_keep = map(str.strip,aopts['GENELIST'].strip().split(","))

		logger.info("Considering %i genes/groups in total" % len(genes_to_keep))

		chrom_to_keep = []
		gene_start = []
		gene_end = []

		# Get gene coordinates from groupfile
		with open(final_groupfile_name) as groupfile:
			for grp_line in groupfile:
				line_s = grp_line.rstrip().split()
				line_gene = line_s[0]

				if line_gene in genes_to_keep:
					parsed_variants = map(parse_epacts,line_s[1:])

					# All of the variants should have the same chromosome. If they don't, something is wrong.
					if len(set([x[0] for x in parsed_variants])) > 1:
						die("Error: group file had variants from multiple chromosomes for the same gene. Line was: %s" % grp_line.rstrip(),logger)

					chrtokeep = parsed_variants[0][0] # Just use the chrom of the first variant. They're all the same.
					chrom_to_keep.append(chrtokeep)

					gene_pos = [x[1] for x in parsed_variants]
					gene_start.append(min(gene_pos) - 1)
					gene_end.append(max(gene_pos) + 1)

		# create bedfile with these position
		tobed = pandas.DataFrame({
			'CHR' : chrom_to_keep,
			'START' : gene_start,
			'END' : gene_end
		})

		# Make sure BED columns in correct order (dictionary above re-orders according to the key)
		tobed = tobed.reindex(columns = ["CHR","START","END"])

		# Remove duplicated gene regions (it's just for tabixing)
		tobed = tobed[~ tobed.duplicated()]

		# Sort (for easier glancing at)
		tobed.sort(["CHR","START"],inplace=True)

		gene_bed = aopts['OUTPREFIX'] + '.genes_to_keep.bed'
		gene_vcf = aopts['OUTPREFIX'] + '.genes_to_keep.vcf'

		tobed.to_csv(gene_bed,sep="\t",index=False,header=False)

		logger.info("Extracting markers from original VCFs needed for single/group tests...")

		if not aopts["SEPCHR"]:
			# Extract the positions from the bedfile and create a separate vcf
			tabixcommand = "tabix -h -B {vcf} {bed} | bgzip -c >| {out}".format(
				vcf = orig_vcf_path,
				bed = gene_bed,
				out = gene_vcf + ".gz"
			)
			logger.debug(tabixcommand)
			run_bash(tabixcommand)

			tabixcommand = 'tabix -p vcf -f ' + gene_vcf + ".gz"
			logger.debug(tabixcommand)
			run_bash(tabixcommand)

		else:
			first = True
			for vcf in match_sepchr_vcfs(orig_vcf_path):
				if first:
					# Grab VCF header. Assume the VCF header for the first chromosome VCF is suitable for the remaining
					# separated VCF files.
					tabixcommand = "tabix -H {vcf} > {out}".format(
						vcf = orig_vcf_path,
						out = gene_vcf
					)
					logger.debug(tabixcommand)
					run_bash(tabixcommand)

					first = False

				# Extract the positions from the bedfile and create a separate vcf
				tabixcommand = "tabix -B {vcf} {bed} >> {out}".format(
					vcf = vcf,
					bed = gene_bed,
					out = gene_vcf
				)
				logger.debug(tabixcommand)
				run_bash(tabixcommand)

			bgzip_cmd = "bgzip -c {vcf} >| {out}".format(
				vcf = gene_vcf,
				out = gene_vcf + ".gz"
			)
			logger.debug(bgzip_cmd)
			run_bash(bgzip_cmd)

			tabixcommand = 'tabix -p vcf -f ' + gene_vcf + ".gz"
			logger.debug(tabixcommand)
			run_bash(tabixcommand)

			try:
				os.remove(gene_vcf)
			except:
				pass

		vcf_for_tests = gene_vcf + ".gz"

		# Do we need a kinship matrix?
		need_kinship_group = sets_overlap(tests,KINSHIP_TESTS)
		need_kinship_single = sets_overlap([aopts["SINGLEMARKERTEST"]],KINSHIP_TESTS)
		kinship_needed = need_kinship_group | need_kinship_single

		# Create a kinship matrix if needed.
		final_kinship_file = None
		if kinship_needed:
			if "KINSHIPFILE" not in aopts:
				logger.info("Creating kinship matrix...")

				# I don't think you ever want to use a few genes to calculate the kinship matrix...
				#vcf_for_kinship = vcf_file_path if aopts["GENELIST"] is not None else final_vcf

				# This is the original VCF file without reduction to a specific set of genes (should estimate from all variants)
				vcf_for_kinship = orig_vcf_path
				final_kinship_file = aopts['OUTPREFIX'] + '.kinf'

				make_kin_cmd = "{epacts} make-kin -vcf {vcf} {sepchr} -ped {ped} -out {out} {filters} -run 1"
				make_kin_cmd = make_kin_cmd.format(
					epacts = epacts,
					vcf = vcf_for_kinship,
					ped = aopts['OUTPREFIX'] + '.pheno.ped',
					out = final_kinship_file,
					sepchr = "-sepchr" if aopts["SEPCHR"] else "",
					filters = EPACTS_KIN_FILTER_ARGS
				)

				run_bash(make_kin_cmd)
				logger.debug(make_kin_cmd)

			else:
				final_kinship_file = aopts["KINSHIPFILE"]

		# Run single marker epacts test
		logger.info("Running single variant EPACTS test...")

		epacts_cmd = "{epacts_loc} single --vcf {vcf} -ped {ped} -pheno {pheno} {covar_cmd} -test {test} {kinship_cmd} " \
								 "{marker_min_maf} {marker_max_maf} {marker_min_mac} -out {out} {field} -run {njobs} {mosix} " \
								 "{remlf}"

		epacts_cmd = epacts_cmd.format(
			epacts_loc = epacts,
			vcf = vcf_for_tests,
			ped = aopts['OUTPREFIX'] + '.pheno.ped',
			pheno = phenotype,
			covar_cmd = covar_cmd,
			test = aopts['SINGLEMARKERTEST'],
			kinship_cmd = "-kin %s" % final_kinship_file if need_kinship_single else "",
			marker_min_maf = marker_min_maf,
			marker_max_maf = marker_max_maf,
			marker_min_mac = marker_min_mac,
			field = field,
			out = aopts['OUTPREFIX'] + '.singlemarker',
			njobs = aopts["NJOBS"],
			mosix = aopts.get("MOSIX",""),
			remlf = "-remlf %s" % aopts["REMLFILE"] if "REMLFILE" in aopts else ""
		)

		logger.debug(epacts_cmd)
		run_bash(epacts_cmd)

		# Read in single marker results
		single_marker_results = pandas.read_table(aopts['OUTPREFIX'] + '.singlemarker.epacts.gz',compression="gzip")
		single_marker_results["MARKER_ID"] = single_marker_results["MARKER_ID"].map(epacts_no_extra)
		single_marker_results = single_marker_results["MARKER_ID PVALUE BETA NS MAF".split()]
		single_marker_results['MAC'] = single_marker_results['NS'] * 2 * single_marker_results['MAF']
		single_marker_results.index = single_marker_results["MARKER_ID"]
		
		# Get variants that pass the MAF/MAC filter
		pass_snps = single_marker_results[single_marker_results['MAF'] < float(aopts['FILTERMAF'])]
		pass_snps = pass_snps[pass_snps['MAF'] >= aopts['MINMAF']]
		pass_snps = pass_snps[pass_snps['MAC'] >= aopts['MARKERMINMAC']]
		macs = pass_snps[['MARKER_ID','MAC']].drop_duplicates()
		pass_snps = list(pass_snps['MARKER_ID'])

		# Get all markers belonging to a gene
		groupfile = open(final_groupfile_name)
		out_groupfile = aopts["OUTPREFIX"] + ".groupfile.filtered.txt"
		newgroupfile = open(out_groupfile, 'w')
		marker_list_for_genes = dict()
		macs_of_genes = dict()
		genes_passing_filters = set()
		marker_names = []
		
		logger.info("Filtering group file...")

		# read in groupfile, and get the list of variants belonging to each gene. Store this in markerlistforgenes.
		# Also get the genemac and the numvars for each gene and get the genes_passing_filters
		for line in groupfile:
			line = line.rstrip()
			lsplit = line.split()
			originaltemp = line.split()
			genename = lsplit[0]
			lsplit = lsplit[1:]
			originaltemp = originaltemp[1:]
			numvars = 0
			mac = 0
			towrite = genename
			passed = False

			# iterating over all markers in gene, find if gene passes filters
			for ind in range(0,len(lsplit)):
				original = lsplit[ind]
				if lsplit[ind] in pass_snps:
					passed = True
					towrite = towrite + "\t" + original
					numvars += 1
					mac += float(macs[macs['MARKER_ID'] == lsplit[ind]]['MAC'])
			
			if (mac >= aopts['GENEMINMAC']) and (numvars >= aopts['MINVARS']):
				genes_passing_filters.add(genename)
				for ind in range(0,len(lsplit)):
					if lsplit[ind] in pass_snps:
						marker_names.append(originaltemp[ind])
			
			if passed:
				newgroupfile.write(towrite + "\n")

			marker_list_for_genes[genename] = lsplit
			macs_of_genes[genename] = mac
			
		groupfile.close()
		newgroupfile.close()

		# Write out the genes passing filters
		with open(aopts['OUTPREFIX'] + '.genes_passing_filters.txt','w') as out:
			print >> out, "GENE"
			for gene in genes_passing_filters:
				print >> out, gene

		# Run the group based tests.
		for grp_test in tests:
			atype, avalue = grp_test.split("=")

			if 'skat' in avalue.lower() and aopts["SKATO"]:
				skato_cmd = "--skat-o"
			else:
				skato_cmd = ""

			epacts_cmd = "{epacts} {type} -test {test} {skato} -vcf {vcf} -pheno {pheno} -ped {ped} -groupf {groupfile} -out {out} " \
									 "{min_maf} {min_mac} {covar} {user_cmd} {kin_cmd} -run {njobs} {mosix} {remlf}"

			run_cmd = epacts_cmd.format(
				epacts = epacts,
				type = atype,
				test = avalue,
				vcf = vcf_for_tests,
				pheno = phenotype,
				ped = aopts['OUTPREFIX'] + '.pheno.ped',
				groupfile = out_groupfile,
				out = aopts['OUTPREFIX'] + '.' + avalue,
				min_maf = min_maf,
				min_mac = marker_min_mac,
				covar = covar_cmd,
				user_cmd = user_cmd,
				kin_cmd = "-kin %s" % final_kinship_file if need_kinship_group else "",
				njobs = aopts["NJOBS"],
				mosix = aopts.get("MOSIX",""),
				skato = skato_cmd,
				remlf = "-remlf %s" % aopts["REMLFILE"] if "REMLFILE" in aopts else ""
			)

			logger.debug(run_cmd)
			run_bash(run_cmd)

		# Read in group test results
		all_grp_test_results = None
		for grp_test in tests:
			grp_results = pandas.read_table(aopts['OUTPREFIX'] + '.' + grp_test.split('=')[1] + '.epacts',sep="\t",header=0)	# reading in output file for test

			if 'NUM_PASS_VARS' in grp_results.columns:
				grp_results.rename(columns={'NUM_PASS_VARS':'PASS_MARKERS'}, inplace=True)

			if 'BEGIN' not in grp_results.columns:
				grp_results.rename(columns={'BEG':'BEGIN'}, inplace=True)
			
			grp_results = grp_results.sort(['#CHROM','BEGIN'])

			if grp_test == tests[0]:
				all_grp_test_results = grp_results
			else:
				all_grp_test_results = all_grp_test_results.append(grp_results)
		
		# Output dataframe contains concatenated results from all EPACTS tests run
		# Note in this case, since it's gene based tests, the MARKER_ID is actually the "group" or the gene tested
		# in the form of chr:start-end_gene
		all_grp_test_results = all_grp_test_results.sort(["#CHROM",'BEGIN'],inplace=False)
		all_grp_test_results = all_grp_test_results[all_grp_test_results.PVALUE <= float(aopts['PVALUETHRESHOLD'])]
		genes = all_grp_test_results['MARKER_ID'].drop_duplicates()

		# if output is empty, then ??
		if len(all_grp_test_results) == 0:
			logger.warning("EPACTS groupwise test returned no significant results.")
			continue
		
		marker_names = pandas.DataFrame({ "MARKER_ID" : marker_names })

		marker_names["CHROM"] = marker_names.iloc[:,0].map(lambda x: x.split(":")[0])
		marker_names["POS"] = marker_names.iloc[:,0].map(lambda x: x.split("_")[0].split(":")[1]).astype("int")
		marker_names["START"] = marker_names["POS"] - 1
		marker_names["END"] = marker_names["POS"] + 1
		marker_names.sort(["CHROM","POS"],inplace=True)

		sig_genes_bed = aopts['OUTPREFIX'] + '.variants_from_sig_genes.txt'

		marker_names_bed = marker_names["CHROM START END".split()]
		marker_names_bed.to_csv(sig_genes_bed,sep="\t",index=False,index_label=False,header=False)

		# Create a VCF file with only the variants from within significant genes
		sig_genes_vcf = aopts['OUTPREFIX'] + '.variants_from_sig_genes.recode.vcf'

		try:
			os.remove(sig_genes_vcf)
		except:
			pass

		# Write VCF header
		with open(sig_genes_vcf,"w") as out:
			orig_header = "\t".join(get_vcf_header(orig_vcf_path))
			print >> out, orig_header

		# Extract variants within significant genes only and write them to a VCF
		tabix_cmd = "tabix -B {vcf} {bed} >> {outvcf}".format(
			vcf = vcf_for_tests,
			bed = sig_genes_bed,
			outvcf = sig_genes_vcf
		)

		logger.debug(tabix_cmd)
		run_bash(tabix_cmd)

		df_sig_gene_vars = pandas.read_table(sig_genes_vcf,sep="\t")
		df_sig_gene_vars.rename(columns = {"#CHROM" : "CHROM"},inplace=True)

		# Remove samples that were filtered out in the PED file
		drop_samples = set(df_sig_gene_vars.columns[9:]).difference(keep_samples)
		df_sig_gene_vars.drop(drop_samples,axis=1,inplace=True)

		# Now for each significant gene, get the genotypes for each sample.
		# We need to melt the genotypes first (into rows for each marker * sample combination), e.g.:
		#  MARKER_ID        IND        GENOTYPE
		#  18:47091686_G/A  MBDT-901   0/0
		#  18:47095862_C/T  MBDT-901   0/0
		#  18:47101838_G/A  MBDT-901   0/0
		#  18:47109939_G/A  MBDT-901   0/0
		#  18:47109955_A/G  MBDT-901   0/0

		for_melt = df_sig_gene_vars.iloc[:,9:]																													# Extract only the sample columns
		for_melt.insert(0,"MARKER_ID",epacts_id_from_df(df_sig_gene_vars))															# Insert MARKER_ID column to melt on
		vcf_sig_gene_melted = pandas.melt(for_melt,"MARKER_ID",var_name="IND",value_name="GENOTYPE")		# Melt the data (see above for example)
		vcf_sig_gene_melted.GENOTYPE = vcf_sig_gene_melted.GENOTYPE.str.replace(':.*','')							  # Drop remaining fields after GT
		vcf_sig_gene_melted = vcf_sig_gene_melted[~ vcf_sig_gene_melted.GENOTYPE.str.contains("\.")]		# Drop missing genos

		# This part just collects per gene the single variant results along with their genotypes per individual into
		# one large melted data frame.
		allvariants = None
		for gene in genes:
			genename = gene.split('_')[1]
			
			if genename not in genes_passing_filters:
				continue

			# Pull out single variant results for this specific gene's variants
			gene_variants = marker_list_for_genes[genename]
			gene_sv_results = single_marker_results[single_marker_results.MARKER_ID.isin(gene_variants)]

			# Merge them together.
			gene_allvar = pandas.merge(gene_sv_results,vcf_sig_gene_melted,on="MARKER_ID",how="left")

			# Include the group/gene
			gene_allvar["GROUP"] = gene

			# Re-arrange columns
			gene_allvar = gene_allvar.reindex(columns = "GROUP MARKER_ID PVALUE GENOTYPE IND BETA MAF MAC".split())

			if allvariants is None:
				allvariants = gene_allvar
			else:
				allvariants = pandas.concat([allvariants,gene_allvar])

		# Load annotation file.
		# In this case, we have to load the entire annotation file, because the user may have specified a group file.
		# If they did, and they accidentally specified an annotation filter, the annotation results would then be filtered,
		# and cause some of the variants in the group file to be missing annotations.
		# If that doesn't make sense, just trust me, load the entire thing.
		if "ANNOTCOLUMNS" in aopts and "ANNOTFILE" in aopts:
			annot_extra_cols = aopts["ANNOTCOLUMNS"].split(",")

			# Load in the annotation file.
			annofile = load_annotation(
				aopts["ANNOTFILE"],
				aopts['ANNOT_CHROM_COL'],
				aopts['ANNOT_POS_COL'],
				aopts['ANNOT_REF_COL'],
				aopts['ANNOT_ALT_COL'],
				aopts['ANNOT_GENE_COL'],
				annot_extra_cols,
				None
			)

			keep_cols = ["EPACTS"] + annot_extra_cols
			annofile = annofile[keep_cols]

			merged = pandas.merge(allvariants,annofile,left_on='MARKER_ID',right_on='EPACTS',how="left")
			del merged['EPACTS']

			allvariants = merged

		allvariants.to_csv(aopts['OUTPREFIX'] + ".allvariants.txt",sep="\t",index=False,index_label=False)

		# Write out a summary of options given.
		with open(aopts['OUTPREFIX'] + ".summary.txt","w") as summary_file:
			print >> summary_file, "DATE\t%s" % time.strftime("%Y-%m-%d %H:%M:%S")

			for k, v in aopts.iteritems():
				if isinstance(v,set) or isinstance(v,list):
					v = ",".join(v)
				else:
					v = str(v)

				print >> summary_file, "%s\t%s" % (k,v)

		# Create PDFs of plots.
		logger.info("Creating output plots for each significant gene...")
		run_plots(aopts,phenotype,covariates,tests_to_write,aopts["PVALUETHRESHOLD"])

if __name__ == "__main__":
	main()
