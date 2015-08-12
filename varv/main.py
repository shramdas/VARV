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

import os, sys, re, getopt, subprocess, gzip, numpy, time
import operator, optparse, logging, pprint, signal, shlex
import pandas
import pandas.computation.ops
import varv
from copy import deepcopy
from varv.util import *
from distutils.dir_util import mkpath
from itertools import repeat, chain
import functools as ft

PROG_NAME = "VARV"
PROG_VERSION = "0.1.2"
PROG_DATE = "03/24/2015"
PROG_AUTHORS = [
	"Shweta Ramdas (sramdas@umich.edu)",
	"Ryan Welch (welchr@umich.edu)"
]
PROG_URL = "https://github.com/shramdas/varv"

pandas.set_option('chained_assignment',None)
pandas.set_option('display.max_colwidth',120)

_VARV_DEBUG = False

EPACTS_KIN_FILTER_ARGS = "-min-maf 0.01 -min-callrate 0.95"
EPACTS_MAKE_GROUP_ARGS = "-nonsyn"
KINSHIP_TESTS = "q.emmax emmax emmaxVT emmaxCMC emmaxSKAT mmskat".split()

def print_program_header():
	prog_string = "%s - %s (%s)" % (PROG_NAME,PROG_VERSION,PROG_DATE)
	max_length = max(map(len,[prog_string,PROG_URL]) + map(len,PROG_AUTHORS))
	max_length += 4; # buffer

	# Top line
	def print_table_line(width=max_length):
		print "".join(chain("+",repeat('-',max_length),"+"))

	print_table_line()

	# Rows
	print "|%s|" % str.center(prog_string,max_length)
	print_table_line()
	for auth in PROG_AUTHORS:
		print "|%s|" % str.center(auth,max_length)
	print_table_line()
	print "|%s|" % str.center(PROG_URL,max_length)
	print_table_line()
	print ""

def get_defaults():
	"""
	Get a dictionary of default values for program options.
	"""

	defaults = dict()

	defaults['INPUTDIR'] = None
	defaults['VCFFILE'] = None
	defaults['FIELD'] = None
	defaults['COVARIATES'] = []
	defaults['PHENOTYPE'] = None
	defaults['MODEL'] = None
	defaults['REMLFILE'] = None
	defaults['GROUPFILE'] = None
	defaults['KINSHIPFILE'] = None
	defaults['ANNOTFILE'] = None
	defaults['PEDFILE'] = None
	defaults['PEDCOLUMNS'] = []
	defaults['PVALUETHRESHOLD'] = 0.001
	defaults['SEPCHR'] = False
	defaults['SINGLEMARKERTEST'] = 'q.linear'
	defaults['GENEMINMAC'] = 5
	defaults['VERBOSE'] = False
	defaults['GENELIST'] = None
	defaults['ANNOTCOLUMNS'] = []
	defaults['ANNOT_GENE_COL'] = "vepGENE"
	defaults['ANNOT_CHROM_COL'] = "CHROM"
	defaults['ANNOT_POS_COL'] = "POS"
	defaults['ANNOT_REF_COL'] = "REF"
	defaults['ANNOT_ALT_COL'] = "ALT"
	defaults["FILTERANNOT"] = None
	defaults['MIN_MAF'] = 0
	defaults['MAX_MAF'] = 1
	defaults['MIN_MAC'] = 1
	defaults['MINVARS'] = 2
	defaults['EPACTS'] = 'epacts'
	defaults['TABIX'] = 'tabix'
	defaults['EPACTSDIR'] = None
	defaults['EPACTSCMD'] = ''
	defaults['NJOBS'] = 1
	defaults['SKATO'] = False

	#defaults['MOSIX'] = '--mosix-nodes="\`/net/fantasia/home/gjun/bin/pick_mainnode\`"'

	return defaults

def run_plots(options,out_prefix,plot_prefix,phenotype):
	log_key = gen_analysis_key(options)
	info = logging.getLogger(log_key).info
	warn = logging.getLogger(log_key).warning
	debug = logging.getLogger(log_key).debug

	b_verbose = _VARV_DEBUG or options["VERBOSE"]

	# These files must have been previously created by the pipeline.
	required_for_plot = map(lambda x: out_prefix + x,[".single_variant_combined.txt"]);

	if all(map(os.path.isfile,required_for_plot)):
		plot_loc = os.path.join(varv.__path__[0],"R/plots.R");

		cmd = "R --vanilla --slave --args {phenotype} {out_prefix} {plot_prefix} < {script}";
		run_cmd = cmd.format(
			out_prefix = out_prefix,
			plot_prefix = plot_prefix,
			phenotype = phenotype,
			script = plot_loc
		)

		run_bash(run_cmd,verbose=b_verbose)

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
	for k in "MODEL OUTPREFIX TEST EPACTS VCFFILE PEDFILE".split():
		if k not in adict:
			raise ValueError, "Error: option '%s' must be set in your config file" % k

		v = adict[k]

		# Value should not be none for required options
		if v is None:
			raise ValueError, "Error: option '%s' must have a sensible value, got: %s" % (k,str(v))

		# If a required option is a list, it should have at least 1 entry
		if hasattr(v,"__iter__") and len(v) == 0:
			raise ValueError, "Error: option '%s' must have at least one value, got: %s" % (k,str(v))

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

				current["INPUTDIR"] = value

				continue

			if key == "VCFFILE":
				if '/' in value or '\\' in value:
					check_file_exists(value)
					current["VCFFILE"] = value
				else:
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

			if key == "ANNOTCOLUMNS":
				value = map(str.strip,value.strip().split(","))

			if key in "MIN_MAF MAX_MAF".split():
				value = float(value)
				if value < 0 or value > 1:
					die("Error: invalid %s in config file: %s" % (key,str(value)))

			if key == "MIN_MAC":
				value = int(value)
				if value < 0:
					die("Error: invalid %s in config file: %s" % (key,str(value)))

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

			if key == "TABIX":
				tb_path, tb_exe = os.path.split(value)
				if tb_path == '':
					if which(tb_exe) is None:
						die("Error: could not find tabix on path, must specify full path to tabix executable using 'TABIX' in config")
				else:
					if not os.path.isfile(value):
						die("Error: could not find tabix as specified by: %s, must specify full path or simply 'tabix' if it is on your PATH" % value)

			if key == "SEPCHR":
				value = True

			if key == "VERBOSE":
				if value in ("OFF","0","FALSE"):
					value = False
				else:
					value = True

			if key == "SKATO":
				value = True

			if key == "PVALUETHRESHOLD":
				value = float(value)
				if value > 1 or value < 0:
					die("Error: invalid p-value threshold specified in config file: " % str(value))

			if key == "PEDCOLUMNS":
				value = value.split(",")

			# This is apparently a "reset" back to default
			if value in ("NA","NaN","None",".","RESET"):
				value = defaults[key]

			current[key] = value

	return analyses

def gen_analysis_key(an):
	import hashlib

	return hashlib.md5("".join([str(x) for x in an.itervalues()])).hexdigest()

def main(arg_string=None):
	print_program_header()

	op = optparse.OptionParser()
	op.add_option("--plotonly",help="Only create the plot PDF, assuming you've already run the pipeline once.",default=False,action="store_true")

	if arg_string is None:
		opts, args = op.parse_args()
	else:
		opts, args = op.parse_args(shlex.split(arg_string))

	if len(args) == 0:
		die("Error: must specify config file, e.g. varv config.cfg")

	if not os.path.isfile(args[0]):
		die("Error: config file does not exist: %s" % args[0])

	# Read in the config file into a list of dictionaries, each one specifying a particular analysis to run
	analyses = read_config(args[0])

	# Loop over each analysis specified in the config file. An analysis is just a set of instructions between PROCESS
	# commands.
	for aopts in analyses:
		tests = aopts["TEST"]
		b_verbose = aopts["VERBOSE"] or _VARV_DEBUG

		# Freeze the verbose argument of run_bash to be dependent on whether
		# we're in debug or verbose mode.
		run_bash_vlevel = ft.partial(run_bash,verbose=b_verbose)

		if len(tests) == 0:
			die("Error: must specify at least 1 test for each PROCESS.")

		outprefix = aopts["OUTPREFIX"]
		if "NETCWD" in outprefix:
			outprefix = re.sub("{{\s*NETCWD\s*}}",net_cwd(),outprefix)

		outprefix_dir = os.path.split(outprefix)[0]
		if not os.path.isdir(outprefix_dir):
			print "Creating directory: %s" % outprefix_dir
			mkpath(outprefix_dir)

		print "Output directory: %s" % outprefix_dir

		# A unique key for this particular analysis. Just makes it easier to setup a logger.
		analysis_key = gen_analysis_key(aopts)

		# Create a logger to log to both STDOUT and a log file
		logger = logging.getLogger(analysis_key)

		# Setup logging to file
		ffmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
		fhandler = logging.FileHandler(outprefix + ".log")
		fhandler.setFormatter(ffmt)
		logger.addHandler(fhandler)

		# Setup logging to STDOUT
		logger.addHandler(logging.StreamHandler(sys.stdout))

		# If we're in debug mode, or they want verbose output, print everything.
		if b_verbose:
			logger.setLevel(logging.DEBUG)
		else:
			logger.setLevel(logging.INFO)

		logger.info("Running analysis: ")
		df_opts = pandas.DataFrame([(k,v) for k,v in aopts.iteritems()],columns=["Option","Value"])
		df_opts.sort("Option",inplace=True)
		logger.info(df_opts.to_string(index=False))
		logger.info("")

		orig_vcf_path = aopts["VCFFILE"]
		epacts = aopts["EPACTS"]
		tabix = aopts["TABIX"]
		pedcolumns = aopts["PEDCOLUMNS"]
		user_cmd = aopts.get("EPACTSCMD","")
		field = '-field %s' % aopts["FIELD"] if aopts.get("FIELD") is not None else ' '
		min_maf = '-min-maf ' + str(aopts['MIN_MAF'])
		max_maf = '-max-maf ' + str(aopts['MAX_MAF'])
		min_mac = '-min-mac ' + str(aopts['MIN_MAC'])
		covariates = aopts["COVARIATES"]
		phenotype = aopts["PHENOTYPE"]
		group_pval_threshold = aopts["PVALUETHRESHOLD"]

		tests_to_write = ",".join([t.split("=")[1] for t in tests])

		if opts.plotonly:
			#logger.info("Creating plots...")
			#run_plots(aopts,phenotype,covariates,tests_to_write,aopts["PVALUETHRESHOLD"])

			raise NotImplementedError, "--plotonly temporarily removed while new plotting code underway"

		# First step: Index VCF if tabix index file does not already exist
		# if sepchr is on, then tabix all other chromosomes too
		if aopts["SEPCHR"]:
			sep_vcfs = match_sepchr_vcfs(orig_vcf_path)
			for f in sep_vcfs:
				if os.path.isfile(f + ".tbi"):
					continue

				tabix_cmd = "%s -p vcf %s" % (tabix,f)
				logger.debug(tabix_cmd)
				run_bash_vlevel(tabix_cmd)

		elif not os.path.isfile(orig_vcf_path + '.tbi'):
			logger.info("Indexing VCF...")

			tabix_command = '%s -p vcf -f %s' % (tabix,orig_vcf_path)
			logger.debug(tabix_command)
			run_bash_vlevel(tabix_command)

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

		pedfile.iloc[:,1].to_csv(outprefix + ".samples_to_keep.txt", sep="\t",index=False, na_rep='NA')
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

		pedfile.to_csv(outprefix + '.pheno.ped', sep="\t",index=False, na_rep='NA')

		# Create groupfile if one is not already given by the user
		if aopts.get("GROUPFILE") is None:
			logger.info("Creating groupfile...")

			final_groupfile_name = outprefix + '.groupfile.txt'

			# No annotation file provided by the user, so we have to use epacts to create the groupfile
			if aopts.get("ANNOTFILE") is None:
				logger.info("No annotation file provided, using EPACTS to create group file")

				if not aopts["SEPCHR"]:
					annot_cmd = "{epacts} anno -in {vcf} -out {out}".format(
						epacts = epacts,
						vcf = orig_vcf_path,
						out = outprefix + '.anno.vcf.gz'
					)

					logger.debug(annot_cmd)
					run_bash_vlevel(annot_cmd)

					annot_cmd = "{epacts} make-group -vcf {vcf} -out {out} --format epacts {args}".format(
						epacts = epacts,
						vcf = outprefix + '.anno.vcf.gz',
						out = final_groupfile_name,
						args = EPACTS_MAKE_GROUP_ARGS
					)
					logger.debug(annot_cmd)
					run_bash_vlevel(annot_cmd)

				# SEPCHR, so we annotate each chromosome separately
				else:
					sep_vcfs = match_sepchr_vcfs(aopts["VCFFILE"])

					try:
						safe_delete_file(final_groupfile_name)
					except:
						pass

					# Loop over each of the VCF files separated by chromosome
					for sep_vcf_fpath in sep_vcfs:
						sep_vcf_name = os.path.split(sep_vcf_fpath)[1]

						annot_cmd = "{epacts} anno -in {vcf} -out {out}".format(
							epacts = epacts,
							vcf = sep_vcf_fpath,
							out = outprefix + '.anno.' + sep_vcf_name
						)
						logger.debug(annot_cmd)
						run_bash_vlevel(annot_cmd)

						annot_cmd = "{epacts} make-group -vcf {vcf} -out {out} -format epacts {args}".format(
							epacts = epacts,
							vcf = outprefix + '.anno.' + sep_vcf_name,
							out = outprefix + '.group.' + sep_vcf_name,
							args = EPACTS_MAKE_GROUP_ARGS
						)
						logger.debug(annot_cmd)
						run_bash_vlevel(annot_cmd)

						annot_cmd = "cat {0} >> {1}".format(
							outprefix + '.group.' + sep_vcf_name,
							final_groupfile_name
						)
						logger.debug(annot_cmd)
						run_bash_vlevel(annot_cmd)

			# They specified an annotation file, so we should use that instead of EPACTS' anno.
			else:
				# Load in the annotation file.
				# If FILTERANNOT was specified, it will be used to filter the data frame.
				annofile = load_annotation(
					aopts["ANNOTFILE"],
					aopts['ANNOT_CHROM_COL'],
					aopts['ANNOT_POS_COL'],
					aopts['ANNOT_REF_COL'],
					aopts['ANNOT_ALT_COL'],
					aopts['ANNOT_GENE_COL'],
					aopts['ANNOTCOLUMNS'],
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

		gene_bed = outprefix + '.genes_to_keep.bed'
		gene_vcf = outprefix + '.genes_to_keep.vcf'

		tobed.to_csv(gene_bed,sep="\t",index=False,header=False)
		run_bash_vlevel('sort -g -k1 -k2 ' + gene_bed + ' > ' + gene_bed + '.sorted')
		run_bash_vlevel('cp ' + gene_bed + '.sorted  ' + gene_bed)

		logger.info("Extracting markers from original VCFs needed for single/group tests...")

		if not aopts["SEPCHR"]:
			# Extract the positions from the bedfile and create a separate vcf
			tabixcommand = "{tabix} -h -R {bed} {vcf} | bgzip -c >| {out}".format(
				tabix = tabix,
				vcf = orig_vcf_path,
				bed = gene_bed,
				out = gene_vcf + ".gz"
			)
			logger.debug(tabixcommand)
			run_bash_vlevel(tabixcommand)
	
			#sort vcf
			run_bash_vlevel("zcat " + gene_vcf + '.gz | grep  "#" > ' + outprefix + 'temp')
			run_bash_vlevel("zcat " + gene_vcf + '.gz | grep -v  "#" | sort -g -k1 -k2 >> ' + outprefix + 'temp')
			
			run_bash_vlevel('mv ' + outprefix + 'temp ' + gene_vcf)
			run_bash_vlevel('rm ' + gene_vcf + '.gz')

			
			run_bash_vlevel('bgzip -c ' + gene_vcf + " >| " + gene_vcf + ".gz")

			tabixcommand = "%s -p vcf -f %s" % (tabix,gene_vcf + ".gz")
			logger.debug(tabixcommand)
			run_bash_vlevel(tabixcommand)

		else:
			first = True
			for vcf in match_sepchr_vcfs(orig_vcf_path):
				if first:
					# Grab VCF header. Assume the VCF header for the first chromosome VCF is suitable for the remaining
					# separated VCF files.
					tabixcommand = "{tabix} -H {vcf} > {out}".format(
						tabix = tabix,
						vcf = orig_vcf_path,
						out = gene_vcf
					)
					logger.debug(tabixcommand)
					run_bash_vlevel(tabixcommand)

					first = False

				# Extract the positions from the bedfile and create a separate vcf
				tabixcommand = "{tabix} -R {bed} {vcf} >> {out}".format(
					tabix = tabix,
					vcf = vcf,
					bed = gene_bed,
					out = gene_vcf
				)
				logger.debug(tabixcommand)
				run_bash_vlevel(tabixcommand)

			bgzip_cmd = "bgzip -c {vcf} >| {out}".format(
				vcf = gene_vcf,
				out = gene_vcf + ".gz"
			)
			logger.debug(bgzip_cmd)
			run_bash_vlevel(bgzip_cmd)

			tabixcommand = "%s -p vcf -f %s" % (tabix,gene_vcf + ".gz")
			logger.debug(tabixcommand)
			run_bash_vlevel(tabixcommand)

			try:
				safe_delete_file(gene_vcf)
			except:
				pass

		vcf_for_tests = gene_vcf + ".gz"
		
		# Do we need to drop samples from the VCF as well?
		vcf_samples = set(vcf_get_header(vcf_for_tests)[9:])
		
		keep_samples = keep_samples[keep_samples.isin(vcf_samples)]
		keep_samples = pandas.Series(keep_samples)
		
		if len(vcf_samples) > 0:
			logger.info("Reducing samples in VCF to match PED file...")

			# TODO: fix bgzip pathing here (and everywhere else)
			vcf_filter_samples(vcf_for_tests,keep_samples,tabix)

		# Do we need a kinship matrix?
		grouptests = []
		for i in aopts['TEST']:
			grouptests.append(i.replace('group=','').strip())
		
		grouptests = set(grouptests)
		need_kinship_group = sets_overlap(grouptests,KINSHIP_TESTS)
		need_kinship_single = sets_overlap([aopts["SINGLEMARKERTEST"]],KINSHIP_TESTS)
		kinship_needed = need_kinship_group | need_kinship_single

		# Create a kinship matrix if needed.
		final_kinship_file = None
		if kinship_needed:
			if aopts.get("KINSHIPFILE") is None:
				logger.info("Creating kinship matrix...")

				# I don't think you ever want to use a few genes to calculate the kinship matrix...
				#vcf_for_kinship = vcf_file_path if aopts["GENELIST"] is not None else final_vcf

				# This is the original VCF file without reduction to a specific set of genes (should estimate from all variants)
				vcf_for_kinship = orig_vcf_path
				final_kinship_file = outprefix + '.kinf'

				make_kin_cmd = "{epacts} make-kin -vcf {vcf} {sepchr} -ped {ped} -out {out} {filters} -run 1"
				make_kin_cmd = make_kin_cmd.format(
					epacts = epacts,
					vcf = vcf_for_kinship,
					ped = outprefix + '.pheno.ped',
					out = final_kinship_file,
					sepchr = "-sepchr" if aopts["SEPCHR"] else "",
					filters = EPACTS_KIN_FILTER_ARGS
				)

				run_bash_vlevel(make_kin_cmd)
				logger.debug(make_kin_cmd)

			else:
				final_kinship_file = aopts["KINSHIPFILE"]

		# Run single marker epacts test
		logger.info("Running single variant EPACTS...")

		epacts_cmd = "{epacts_loc} single --vcf {vcf} -ped {ped} -pheno {pheno} {covar_cmd} -test {test} {kinship_cmd} " \
								 "{marker_min_maf} {marker_max_maf} {marker_min_mac} -out {out} {field} -run {njobs} {mosix} " \
								 "{remlf}"

		epacts_cmd = epacts_cmd.format(
			epacts_loc = epacts,
			vcf = vcf_for_tests,
			ped = outprefix + '.pheno.ped',
			pheno = phenotype,
			covar_cmd = covar_cmd,
			test = aopts['SINGLEMARKERTEST'],
			kinship_cmd = "-kin %s" % final_kinship_file if need_kinship_single else "",
			marker_min_maf = min_maf,
			marker_max_maf = max_maf,
			marker_min_mac = min_mac,
			field = field,
			out = outprefix + '.singlemarker',
			njobs = aopts["NJOBS"],
			mosix = aopts.get("MOSIX",""),
			remlf = "-remlf %s" % aopts["REMLFILE"] if aopts.get("REMLFILE") is not None else ""
		)

		logger.debug(epacts_cmd)
		run_bash_vlevel(epacts_cmd)

		# Read in single marker results.
		single_marker_results = pandas.read_table(outprefix + '.singlemarker.epacts.gz',compression="gzip")
		single_marker_results["MARKER_ID"] = single_marker_results["MARKER_ID"].map(epacts_no_extra)
		single_marker_results = single_marker_results["MARKER_ID PVALUE BETA NS AC MAF".split()]
		single_marker_results["MAC"] = single_marker_results.apply(lambda x: min(2 * x["NS"] - x["AC"],x["AC"]),axis=1)
		single_marker_results.index = single_marker_results["MARKER_ID"]
		
		# Get variants that pass the MAF/MAC filter
		pass_snps = single_marker_results[single_marker_results['MAF'] < aopts['MAX_MAF']]
		pass_snps = pass_snps[pass_snps['MAF'] >= aopts['MIN_MAF']]
		pass_snps = pass_snps[pass_snps['MAC'] >= aopts['MIN_MAC']]
		macs = pass_snps[['MARKER_ID','MAC']].drop_duplicates()
		pass_snps = list(pass_snps['MARKER_ID'])

		# Get all markers belonging to a gene
		groupfile = open(final_groupfile_name)
		out_groupfile = outprefix + ".groupfile.filtered.txt"
		newgroupfile = open(out_groupfile, 'w')
		marker_list_for_genes = dict()
		genes_passing_filters = set()

		logger.info("Filtering group file...")

		with groupfile, newgroupfile:
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
				else:
					#logger.debug("gene %s failed filters - GENEMINMAC: %f | MINVARS: %i" % (genename,mac,numvars))
					pass

				if passed:
					newgroupfile.write(towrite + "\n")

				marker_list_for_genes[genename] = lsplit

		# Write out the genes passing filters
		with open(outprefix + '.genes_passing_filters.txt','w') as out:
			print >> out, "GENE"
			for gene in genes_passing_filters:
				print >> out, gene

		# Run the group based tests.
		logger.info("Running gene-based tests...")
		for grp_test in tests:
			atype, avalue = grp_test.split("=")

			if 'skat' in avalue.lower() and aopts["SKATO"]:
				skato_cmd = "--skat-o"
			else:
				skato_cmd = ""

			epacts_cmd = "{epacts} {type} -test {test} {skato} -vcf {vcf} -pheno {pheno} -ped {ped} -groupf {groupfile} -out {out} " \
									 "{min_maf} {min_mac} {max_maf} {covar} {user_cmd} {kin_cmd} -run {njobs} {mosix} {remlf}"

			run_cmd = epacts_cmd.format(
				epacts = epacts,
				type = atype,
				test = avalue,
				vcf = vcf_for_tests,
				pheno = phenotype,
				ped = outprefix + '.pheno.ped',
				groupfile = out_groupfile,
				out = outprefix + '.' + avalue,
				min_maf = min_maf,
				min_mac = min_mac,
				max_maf = max_maf,
				covar = covar_cmd,
				user_cmd = user_cmd,
				kin_cmd = "-kin %s" % final_kinship_file if need_kinship_group else "",
				njobs = aopts["NJOBS"],
				mosix = aopts.get("MOSIX",""),
				skato = skato_cmd,
				remlf = "-remlf %s" % aopts["REMLFILE"] if aopts.get("REMLFILE") is not None else ""
			)

			logger.debug(run_cmd)
			run_bash_vlevel(run_cmd)

		logger.info("Starting to create output files and plots...")

		# Read in group test results
		all_grp_test_results = None
		for grp_i, grp_test in enumerate(tests):
			grp_test_name = grp_test.split("=")[1]
			grp_results = pandas.read_table(outprefix + '.' + grp_test_name + '.epacts',sep="\t",header=0)	# reading in output file for test

			if 'NUM_PASS_VARS' in grp_results.columns:
				grp_results.rename(columns={'NUM_PASS_VARS':'PASS_MARKERS'}, inplace=True)

			if 'BEGIN' not in grp_results.columns:
				grp_results.rename(columns={'BEG':'BEGIN'}, inplace=True)
			
			grp_results = grp_results.sort(['#CHROM','BEGIN'])

			# The name of the gene-based test
			grp_results["TEST"] = grp_test_name

			# Extract the name of the gene from the EPACTS MARKER_ID column
			grp_results["GENE"] = grp_results.MARKER_ID.map(lambda x: x.split("_")[-1])

			if grp_i == 0:
				all_grp_test_results = grp_results
			else:
				all_grp_test_results = all_grp_test_results.append(grp_results)
		
		# Output dataframe contains concatenated results from all EPACTS tests run
		# Note in this case, since it's gene based tests, the MARKER_ID is actually the "group" or the gene tested
		# in the form of chr:start-end_gene
		all_grp_test_results = all_grp_test_results.sort(["#CHROM",'BEGIN'],inplace=False)

		# Make a separate df with only the significant gene based test results
		all_sig_grp_test_results = all_grp_test_results[all_grp_test_results.PVALUE <= group_pval_threshold]

		# List of genes that were significant
		sig_groups = all_sig_grp_test_results['MARKER_ID'].drop_duplicates()
		sig_genes = map(lambda x: x.split("_")[1],sig_groups)

		# if output is empty, then ??
		if len(all_sig_grp_test_results) == 0:
			logger.warning("EPACTS groupwise test returned no significant results.")
			continue
	
		sig_gene_markers = reduce(operator.add,map(marker_list_for_genes.get,all_sig_grp_test_results.GENE))
		sig_gene_markers = filter(lambda x: x is not None,sig_gene_markers)
		sig_gene_markers = set(sig_gene_markers)

		for_sig_marker_bed = pandas.DataFrame({ "MARKER_ID" : list(sig_gene_markers) })
		for_sig_marker_bed["CHROM"] = for_sig_marker_bed.iloc[:,0].map(lambda x: x.split(":")[0])
		for_sig_marker_bed["POS"] = for_sig_marker_bed.iloc[:,0].map(lambda x: x.split("_")[0].split(":")[1]).astype("int")
		for_sig_marker_bed["START"] = for_sig_marker_bed["POS"] - 1
		for_sig_marker_bed["END"] = for_sig_marker_bed["POS"] + 1
		for_sig_marker_bed.sort(["CHROM","POS"],inplace=True)

		sig_genes_bed = outprefix + '.variants_from_sig_genes.txt'
		marker_names_bed = for_sig_marker_bed["CHROM START END".split()]
		marker_names_bed.to_csv(sig_genes_bed,sep="\t",index=False,index_label=False,header=False)

		# Create a VCF file with only the variants from within significant genes
		sig_genes_vcf = outprefix + '.variants_from_sig_genes.recode.vcf'

		try:
			safe_delete_file(sig_genes_vcf)
		except:
			pass

		# Extract variants within significant genes only and write them to a VCF
		tabix_cmd = "{tabix} -h -R {bed} {vcf} >> {outvcf}".format(
			tabix = tabix,
			vcf = vcf_for_tests,
			bed = sig_genes_bed,
			outvcf = sig_genes_vcf
		)

		logger.debug(tabix_cmd)
		run_bash_vlevel(tabix_cmd)

		df_sig_gene_vars = vcf_pandas_load(sig_genes_vcf)
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

		def merge_group_and_single_results(group_results,single_results,markers_per_gene,passing_genes):
			"""
			This function collects per gene the single variant results along with the gene level test results
			into one final data frame, like this:

				 GENE      mmskat.P          VARIANT      SV.P     BETA      MAF  MAC  GENE_FILTERED  GENE_SIG
			0  LCAT  7.242200e-01  16:67976823_C/T  0.724200  0.04040  0.00467   79              1         0
			1  LIPG  6.196300e-11  18:47091686_G/A  0.000003  0.76050  0.00225   38              0         1
			2  LIPG  6.196300e-11  18:47095862_C/T  0.000184  0.43700  0.00450   76              0         1
			3  LIPG  6.196300e-11  18:47101838_G/A  0.791900 -0.26300  0.00006    1              0         1
			4  LIPG  6.196300e-11  18:47109939_G/A  0.957000 -0.01924  0.00047    8              0         1
			5  LIPG  6.196300e-11  18:47109955_A/G  0.000044  0.36990  0.00728  123              0         1
			"""

			pivot_grp_results = group_results["MARKER_ID TEST PVALUE".split()].pivot("MARKER_ID","TEST","PVALUE")
			pivot_grp_results.rename(columns = dict(zip(pivot_grp_results.columns,map(lambda x: x + ".P",pivot_grp_results.columns))),inplace=True)
			group_p_cols = list(pivot_grp_results.columns)
			pivot_grp_results["GENE_SIG"] = pivot_grp_results.apply(lambda x: (x < group_pval_threshold).any(),axis=1).astype("int")

			# This crazy thing just takes a dictionary of gene --> [marker list] pairs and converts it into a unpivoted data frame, like:
			# gene1 rs1
			# gene1 rs2
			# gene1 rs3
			# gene2 rs8
			df_markers_for_genes = pandas.DataFrame.from_records(reduce(operator.add,map(lambda x: zip(repeat(x[0]),x[1]),markers_per_gene.iteritems())))
			df_markers_for_genes.columns = "GENE VARIANT".split()
			df_markers_for_genes["GENE_FILTERED"] = (~df_markers_for_genes["GENE"].isin(passing_genes)).astype("int")

			gsv1 = group_results["MARKER_ID GENE".split()]
			gsv1.rename(columns = {"MARKER_ID" : "GROUP"},inplace=True)

			gsv2 = pandas.merge(gsv1,df_markers_for_genes,on="GENE")
			gsv3 = pandas.merge(gsv2,single_results,left_on="VARIANT",right_on="MARKER_ID")
			gsv3.rename(columns = {
				"PVALUE" : "SV.P",
			},inplace=True)
			gsv4 = pandas.merge(gsv3,pivot_grp_results,left_on="GROUP",right_index=True)
			gsv5 = gsv4.drop("GROUP MARKER_ID AC NS".split(),axis=1)
			s5_col_order = ["GENE"] + group_p_cols + "VARIANT SV.P BETA MAF MAC GENE_FILTERED GENE_SIG".split()
			gsv5 = gsv5[s5_col_order]

			try:
				gsv5["MAC"] = gsv5["MAC"].astype("int")
			except:
				pass

			return gsv5

		gene_final = merge_group_and_single_results(
			all_grp_test_results,
			single_marker_results,
			marker_list_for_genes,
			genes_passing_filters
		)

		# Merge annotations in with single variant results.
		# In this case, we have to load the entire annotation file, because the user may have specified a group file.
		# If they did, and they accidentally specified an annotation filter, the annotation results would then be filtered,
		# and cause some of the variants in the group file to be missing annotations.
		# If that doesn't make sense, just trust me, load the entire thing.
		if aopts.get("ANNOTFILE") is not None:
			# Load in the annotation file.
			annofile = load_annotation(
				aopts["ANNOTFILE"],
				aopts['ANNOT_CHROM_COL'],
				aopts['ANNOT_POS_COL'],
				aopts['ANNOT_REF_COL'],
				aopts['ANNOT_ALT_COL'],
				aopts['ANNOT_GENE_COL'],
				aopts['ANNOTCOLUMNS'],
				None
			)

			keep_cols = ["EPACTS"] + aopts["ANNOTCOLUMNS"]
			annofile = annofile[keep_cols]

			gene_final = pandas.merge(gene_final,annofile,left_on='VARIANT',right_on='EPACTS',how="left")
			del gene_final['EPACTS']

		# Gene/single results reduced to genes that are significant, and passed filters.
		gene_final_filtered = gene_final.query("(GENE_SIG == 1) & (GENE_FILTERED != 1)")
		del gene_final_filtered["GENE_FILTERED"]
		del gene_final_filtered["GENE_SIG"]

		def genotype_trait_tables(genes,genotypes,ped,phenotype,markers_per_gene,outprefix):
			pcols = ["IND_ID",phenotype]
			tables = []

			for gene in genes:
				markers = markers_per_gene.get(gene)
				if markers is None:
					continue

				s1 = genotypes.query("MARKER_ID in @markers")
				s2 = s1.iloc[:,1:].transpose()
				s2.columns = s1.MARKER_ID
				s3 = pandas.merge(ped[pcols],s2,right_index=True,left_on="IND_ID")
				s3.sort(phenotype,inplace=True)

				tables.append((gene,s3))

				out = outprefix + ".genotype_trait_table.gene-%s.tab" % gene
				s3.to_csv(out,sep="\t",index=False)

			return tables

		gt_tables = genotype_trait_tables(gene_final_filtered.GENE.unique(),for_melt,pedfile,phenotype,marker_list_for_genes,outprefix)

		# Where should we write the top genes plot?
		html_dir = os.path.join(outprefix + ".plots/")
		mkpath(html_dir)

		# Write the supporting files (HTML and JS) for the plot.
		copy_html_template(html_dir)

		# Write out phenotype data for plotting.
		for_plot_pheno = pedfile[[phenotype,"IND_ID"]]
		for_plot_pheno.columns = "TRAIT IND".split()
		df_to_js(for_plot_pheno,"phenos",os.path.join(html_dir,"plot_phenos.js"),float_format="%0.3g",write_tab=True)

		# Write out single variant and gene test results for plotting.
		# We only want those single variant results for genes passing filters, though.
		df_to_js(gene_final_filtered,"genes",os.path.join(html_dir,"plot_genes.js"),float_format="%0.3g",write_tab=True)

		# Write out single variant results for creating QQ plots and manhattan plots.
		for_qq_manhattan = gene_final["GENE GENE_FILTERED VARIANT SV.P".split()]
		for_qq_manhattan["CHROM"] = for_qq_manhattan.VARIANT.map(lambda x: parse_epacts(x)[0])
		for_qq_manhattan["POS"] = for_qq_manhattan.VARIANT.map(lambda x: parse_epacts(x)[1])
		for_qq_manhattan.to_csv(os.path.join(outprefix + ".single_variant_combined.txt"),sep="\t",index=False,na_rep="NA")

		# Write out variant data for plotting.
		for_plot_variants = add_rare_count(vcf_sig_gene_melted,filter=0)
		for_plot_variants.rename(columns = {
			"MARKER_ID" : "VARIANT"
		},inplace=True)
		for_plot_variants = pandas.merge(for_plot_variants,for_plot_pheno,on="IND",how="left")
		df_to_js(for_plot_variants,"variants",os.path.join(html_dir,"plot_variants.js"),float_format="%0.3g",write_tab=True)

		#if trait is binary, add information about number of cases/controls each variant is found in:
		if len(for_plot_variants['TRAIT'].unique()) <=2 :
			allv = gene_final_filtered['VARIANT'].unique()
			hets_controls = np.array([])
			hets_cases = np.array([])
			homs_controls = np.array([])
			homs_cases = np.array([])

			for v in allv:
				c = for_plot_variants[for_plot_variants['VARIANT'] == v]
				n1 = len(c[(c['GENOTYPE'] == '0/1') & (c['TRAIT'] == 1)])
				n2 = len(c[(c['GENOTYPE'] == '0/1') & (c['TRAIT'] == 2)])
				n3 = len(c[(c['GENOTYPE'] == '1/1') & (c['TRAIT'] == 1)])
				n4 = len(c[(c['GENOTYPE'] == '1/1') & (c['TRAIT'] == 2)])
				hets_controls = np.append(hets_controls, n1)
				hets_cases = np.append(hets_cases,n2)
				homs_controls = np.append(homs_controls,n3)
				homs_cases = np.append(homs_cases,n4)
			
			hets_controls = pandas.Series(hets_controls, index=allv,dtype='int64')
			hets_cases = pandas.Series(hets_cases, index=allv,dtype='int64')
			homs_cases = pandas.Series(homs_cases, index=allv,dtype='int64')
			homs_controls = pandas.Series(homs_controls, index=allv,dtype='int64')

			allnums = pandas.concat([hets_cases, homs_cases, hets_controls, homs_controls],axis=1)
			allnums.to_csv("allnums.txt",header=True)
			allnums['VARIANT'] = allnums.index
			allnums.columns = ['HETS_CASES','HOMS_CASES','HETS_CONTROLS','HOMS_CONTROLS','VARIANT']
			plot_genes_binary =  pandas.merge(gene_final_filtered, allnums, left_on='VARIANT',right_on='VARIANT',how='inner')
			df_to_js(plot_genes_binary,"genes",os.path.join(html_dir,"plot_genes_binary.js"),float_format="%0.3g",write_tab=True)


		# Write out the model formula for the plotting HTML code.
		with open(os.path.join(html_dir,"plot_info.js"),"w") as out:
			print >> out, 'model_formula = "%s"' % aopts["MODEL"]

		# Write out a summary of options given.
		with open(outprefix + ".summary.txt","w") as summary_file:
			print >> summary_file, "DATE\t%s" % time.strftime("%Y-%m-%d %H:%M:%S")

			for k, v in aopts.iteritems():
				if isinstance(v,set) or isinstance(v,list):
					v = ",".join(v)
				else:
					v = str(v)

				print >> summary_file, "%s\t%s" % (k,v)

		# Create PDFs of plots.
		run_plots(aopts,outprefix,html_dir,phenotype)

		# from IPython.core.debugger import Tracer
		# debug = Tracer()
		# debug()

		logger.info("Completed @ %s" % time.strftime("%H:%m:%S %Y-%M-%d"))

		return opts, args

if __name__ == "__main__":
	main()
