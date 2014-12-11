#!/usr/bin/env python
import sys, os, subprocess, gzip, signal, pprint, re, base64, shutil
import numpy as np
import pandas
import pandas.computation
import runepacts
from collections import Counter

def copy_html_template(out_dir):
	"""
	Copies HTML and supporting javascript files to an output directory. These files currently contain the plotting
	code for the new top genes graphic.

	:param out_dir: Directory to write files to.
	:return: None
	"""

	data_dir = os.path.join(runepacts.__path__[0],"data/")

	# Copy the required javascript files to the outpath.
	files = os.listdir(data_dir)
	for f in files:
		shutil.copy2(os.path.join(data_dir,f),os.path.join(out_dir,f))

def recode_genotypes_additive(genotypes):
	"""
	Recode genotypes to be additive, counting towards
	the minor allele.

	Genotypes are assumed to be in the format:
	00, 01, 10, 11

	Genotypes will be recoded 0, 1, 2
	With 0 - common homozygote
			 1 - heterozygote
			 2 - rare homozygote

	Returns:
	return[0] - recoded genotypes
	return[0] - whether coding was flipped (the major allele was 1 in the input genotypes)
	"""

	genotypes = map(lambda x: x.replace("/","").replace("|",""),genotypes)
	c0 = sum(map(lambda x: x.count('0'),genotypes))
	c1 = (2 * len(genotypes)) - c0

	flip = c1 > c0

	trans = {
		'00' : 0 if not flip else 2,
		'01' : 1,
		'10' : 1,
		'11' : 2 if not flip else 0,
		'..' : np.nan
	}

	recoded = np.array(map(trans.get,genotypes),dtype=np.float16)

	if np.isnan(recoded).all():
		maf = np.nan
	else:
		maf = np.nansum(recoded) / float(2 * len(recoded))

	return recoded, flip, maf

def add_rare_count(dframe,marker_col="MARKER_ID",geno_col="GENOTYPE",count_col="RARE_COUNT",filter=None):
	"""
	Adds a rare count column denoting the number of rare alleles. This function assumes the data frame being passed in
	is in long format, namely that there is a row per individual per genotype.

	:param dframe: Dataframe of genotypes, one row per individual per genotype
	:param marker_col: Name of marker ID column
	:param geno_col: Name of genotype column
	:param count_col: Name for the created rare count column within the data frame
	:param filter: Remove variants with this rare_count or fewer (so if filter = 1, then rare_count <= 1 are dropped)
	:return: The data frame with a new rare count column named by count_col
	"""

	results = []
	for index, df in dframe.groupby(marker_col):
		# Find the rarest allele
		rare_al, rare_total_count = Counter("".join(df[geno_col]).replace("/","")).most_common()[-1]

		# Add in a column counting the number of rare alleles per person.
		df[count_col] = df[geno_col].str.count(rare_al)

		results.append(df)

	last = pandas.concat(results)
	if filter is not None:
		last = last[last[count_col] > filter]

	return last

def df_to_js(dframe,fname,js_out,write_tab=True,float_format=None):
	"""
	Creates a javascript file with a function that returns this data frame as a base64 encoded string.
	This helps get around a browser limitation in Chrome that will not allow reading data from local files.
	Instead, we just include these javascript files with the data embedded.

	:param dframe: Pandas data frame
	:param fname: Name of the javascript function that will return the data
	:param js_out: Name of the javascript file to write
	:param write_tab: If enabled, also write the tab-delimited tsv file next to the js file.
	"""

	data = dframe.to_csv(sep="\t",index=False,na_rep="NA")

	template = """

	function load_%s() {
		return atob("%s")
	}

	""" % (fname,base64.encodestring(data).replace("\n",""))

	with open(js_out,"w") as out:
		print >> out, template

	if write_tab:
		dframe.to_csv(
			js_out.replace(".js",".tsv"),
			sep="\t",
			index=False,
			float_format=float_format,
			na_rep="NA"
		)

def match_sepchr_vcfs(vcf_prefix):
	"""
	Given the filepath of a VCF file that has been separated by chromosome, for example:
	/path/to/study_vcf_chr1.vcf.gz

	Return a list of all separated by chromosome VCFs, e.g.:
	/path/to/study_vcf_chr1.vcf.gz
	/path/to/study_vcf_chr2.vcf.gz
	/path/to/study_vcf_chr3.vcf.gz

	:param vcf_prefix: A VCF file containing a chr* in the filename.
	:return: A list of all matching VCF files.
	"""

	vcf_dir, vcf_file = os.path.split(vcf_prefix)

	try:
		vcf_prefix = re.search("(.*)chr.*",vcf_file).groups()[0]
	except:
		die("Couldn't find chr# in chromosome separated VCF file, should be something like: /path/to/vcf_chr1.vcf.gz, was: %s" % vcf_file)

	matches = []
	for f in os.listdir(vcf_dir):
		is_vcf = f.endswith(".vcf") or f.endswith(".vcf.gz")
		if f.startswith(vcf_prefix) and is_vcf:
			matches.append(os.path.join(vcf_dir,f))

	return matches

def epacts_no_extra(eid):
	"""
	Given an EPACTS ID, remove the extraneous information from the end and only return
	an ID of the form: CHR:POS_REF/ALT
	:param eid: EPACTS ID
	:return: EPACTS ID (without trailing _*)
	"""

	return "_".join(eid.split("_",2)[0:2])

def parse_epacts(v):
	"""
	Try to parse an EPACTS ID into components.
	Returns None if it couldn't be parsed.
	Otherwise a tuple of:
	(chrom, pos, ref, alt, extra)
	"""

	split = v.split("_");

	# It should have at least 2 elements (chrom:pos, ref/alt)
	if len(split) < 2:
		return None;

	# Split chrom/pos
	try:
		chrom, pos = split[0].split(":");
	except:
		return None;

	# Position should be numeric
	try:
		pos = long(pos);
	except:
		return None;

	# Split the alleles
	try:
		ref, alt = split[1].split("/");
	except:
		return None;

	return [chrom,pos,ref,alt,"_".join(split[2:])];

def run_bash(cmd,verbose=True):
	proc = subprocess.Popen(cmd,shell=True,executable="/bin/bash",preexec_fn=os.setsid,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	try:
		stdout, stderr = proc.communicate()
	except KeyboardInterrupt:
		# We want to make sure any processes spawned are terminated if the user hits CTRL+C.
		# This kills the entire process group (os.killpg) and takes advantage of the fact that
		# os.setsid was passed as the preexec_fn above (this causes the child process to create
		# a new process group)
		os.killpg(proc.pid, signal.SIGTERM)
		raise

	if verbose:
		if stdout != '':
			print stdout

		if stderr != '':
			print >> sys.stderr, stderr

	if proc.returncode != 0:
		raise Exception, "Command failed: %s, error code was: %s. See above for error message." % (cmd,str(proc.returncode))

	return proc.returncode

def epacts_id_from_df(dframe,chrom="CHROM",pos="POS",ref="REF",alt="ALT"):
	return dframe[chrom].map(str) + ":" + dframe[pos].map(str) + "_" + dframe[ref].map(str) + '/' + dframe[alt].map(str)

def groupfile_get_genes(filepath):
	"""
	Return the unique set of genes given an EPACTS group file.
	:param filepath: Path to EPACTS group file
	:return: Set of genes
	"""

	genes = set()
	with open(filepath) as f:
		for line in f:
			ls = line.strip()

			if ls == '':
				continue

			gene = ls.split("\t")[0]
			genes.add(gene)

	return genes

def vcf_get_header(vcf):
	if vcf.endswith(".gz"):
		f = gzip.open(vcf);
	else:
		f = open(vcf);
		
	header = None;
	with f:
		for line in f:
			if line.startswith("#CHROM"):
				header = line.rstrip().split("\t");
				break;
				
	return header;

def vcf_nskip(filepath):
	f = gzip.open(filepath) if filepath.endswith(".gz") else open(filepath);

	with f:
		i = 0
		for line in f:
			if line.startswith("##"):
				i += 1
			else:
				break

	return i;

def vcf_pandas_load(vcf,*args,**kwargs):
	skip = vcf_nskip(vcf);
	comp = "gzip" if vcf.endswith(".gz") else None;
	df = pandas.read_table(vcf,sep="\t",compression=comp,skiprows=skip,*args,**kwargs);

	return df;

def sets_overlap(s1,s2):
	"""
	Quick function to check if any elements of two sets overlap.
	:param s1: Iterable 1
	:param s2: Iterable 2
	:return: Returns true if any item in iterable 1 is also in iterable 2
	"""

	return len(set(s1).intersection(s2)) > 0

def calculate_maf(inputvcf,outprefix,samplestokeep,sepchr=False):
	header = vcf_get_header(inputvcf);

	colstokeep = []
	samplestokeeplist = list(samplestokeep)
	for columnnum in range(10,len(header)):
		if header[columnnum] in samplestokeeplist:
			colstokeep.append(columnnum)	
	if inputvcf.endswith(".gz"):
		f = gzip.open(inputvcf);
	else:
		f = open(inputvcf);

	rows = []
	n_alleles = 2
	for line in f:
		if line.startswith("#"):
			continue;

		# Split the line on tabs (VCF must be tab-delimited)
		ls = line.split("\t");
		ls[-1] = ls[-1].rstrip();

		# Find which element of the genotype field is the genotype itself		
		gt_ind = 0;
		
		# Grab genotypes for only those samples we care about.
		# gt_ind is the index within the genotype field for the genotype itself
		# Remember sometimes that field can contain dosage, genotype likelihoods, etc.
		# Example: GT:EC:DS 0/0:1.3314141:8
		genos = [ls[i].split(":")[gt_ind] for i in colstokeep];

		# Count alleles per genotype.
		count = Counter();
		map(count.update,genos);
		n_chr = count['0'] + count['1'];
		freq_0 = count['0'] / float(n_chr);
		freq_1 = count['1'] / float(n_chr);
		mac = min(count['0'],count['1']);
		maf = min(freq_0,freq_1);

		# Store information on this variant
		chrom, pos, idv, ref, alt = ls[0:5];

		rows.append([chrom,pos,idv,ref,alt,n_alleles,n_chr,freq_0,freq_1,mac,maf]);

	df = pandas.DataFrame(
		rows, columns=["CHROM","POS","ID","REF","ALT","N_ALLELES","N_CHR","FREQ0","FREQ1","MAC","MAF"]
	);

	df.to_csv(outprefix + ".maffile.txt",sep="\t",index=False,index_label=False)
	return df

def get_header_names(filepath,indices,sep="\t"):
	"""
	Given a file and a list of column numbers (0-indexed), return the column header names.
	:param filepath:
	:param indices:
	:param sep:
	:return: List of header names. Elements are None if the index was not valid (e.g. off the end, negative, etc.)
	"""
	with open(filepath) as f:
		header = f.readline().split("\t")
		header[-1] = header[-1].rstrip()

	names = []
	for i in indices:
		try:
			names.append(header[i])
		except:
			names.append(None)

	return names

# def load_annotation(anno_file,chrom_col="CHROM",pos_col="POS",ref_col="REF",alt_col="ALT",gene_col="vepGENE",chrpos_col=None,annot_cols=None,annot_filter=None,logger=None):
# 	# Grab the annotation file.
# 	if chrpos_col is not None:
# 		usecols = [chrpos_col,ref_col,alt_col,gene_col]
# 	else:
# 		usecols = [chrom_col,pos_col,ref_col,alt_col,gene_col]
#
# 	if annot_cols is not None:
# 		usecols.extend(annot_cols)
#
# 	compr = "gzip" if anno_file.endswith(".gz") else ""
# 	df_anno = pandas.read_table(
# 		anno_file,
# 		sep = "\t",
# 		usecols = usecols,
# 		compression = compr,
# 		dtype = {
# 			pos_col : np.uint32
# 		}
# 	)
#
# 	# Check that we actually got valid data
# 	if df_anno is None or df_anno.shape[0] <= 0:
# 		die("Error: annotation file contained no data: %s" % anno_file,logger)
#
# 	# Helper function to split CHR:POS field (if necessary; usually they should be giving a file with CHR and POS columns)
# 	def split_chrpos(cp):
# 		cp = cp.replace("chrom","").replace("chr","")
# 		chrom, pos = cp.split("_")[0].split(":")
#
# 		return chrom, pos
#
# 	if chrpos_col is not None:
# 		df_anno[chrom_col] = df_anno[chrpos_col].map(lambda x: split_chrpos(x)[0])
# 		df_anno[pos_col] = df_anno[chrpos_col].map(lambda x: split_chrpos(x)[1])
#
# 	# Fix data types
# 	df_anno[pos_col] = df_anno[pos_col].astype("int")
#
# 	# Add in EPACTS ID
# 	def epacts(row):
# 		return "{0}:{1}_{2}/{3}".format(
# 			row[chrom_col],
# 			row[pos_col],
# 			row[ref_col],
# 			row[alt_col]
# 		)
#
# 	df_anno["EPACTS"] = df_anno.apply(epacts,axis=1)
#
# 	if annot_filter is not None:
# 		try:
# 			df_anno = df_anno.query(annot_filter)
# 		except pandas.computation.ops.UndefinedVariableError:
# 			die("Error: annotation file filter specified variable that does not exist, filter was: %s" % annot_filter,logger)
#
# 		if df_anno is None or df_anno.shape[0] == 0:
# 			die("Error: annotation file filter '%s' resulted in an empty data set" % annot_filter,logger)
#
# 	return df_anno

def load_annotation(anno_file,chrom_col="CHROM",pos_col="POS",ref_col="REF",alt_col="ALT",gene_col="vepGENE",annot_cols=None,annot_filter=None,logger=None):
	"""

	:param anno_file:
	:param chrom_col: Column name for chromosome
	:param pos_col: Column name for position
	:param ref_col: Column name for reference allele
	:param alt_col: Column name for alternate allele
	:param gene_col: Column name for gene (usually vepGENE)
	:param annot_cols: Extra annotation columns to keep along
	:param annot_filter: Optional filter on the data frame. Specify using simple comparisons, e.g. "(a < b) & (c == 1)"
	:param logger: Optional logger object
	:return:
	"""

	if logger is not None:
		logger.debug("load_annotation called: \n" + pprint.pformat(locals()))

	usecols = [chrom_col,pos_col,ref_col,alt_col,gene_col]

	if annot_cols is not None and len(annot_cols) > 0:
		usecols.extend(annot_cols)

	compr = "gzip" if anno_file.endswith(".gz") else ""
	df_anno = pandas.read_table(
		anno_file,
		sep = "\t",
		usecols = usecols,
		compression = compr,
		dtype = {
			pos_col : np.uint32
		}
	)

	# Check that we actually got valid data
	if df_anno is None or df_anno.shape[0] <= 0:
		die("Error: annotation file contained no data: %s" % anno_file,logger)

	# Fix data types
	df_anno[pos_col] = df_anno[pos_col].astype("int")

	# Add in EPACTS ID
	def epacts(row):
		return "{0}:{1}_{2}/{3}".format(
			row[chrom_col],
			row[pos_col],
			row[ref_col],
			row[alt_col]
		)

	df_anno["EPACTS"] = df_anno.apply(epacts,axis=1)

	if annot_filter is not None:
		try:
			df_anno = df_anno.query(annot_filter)
		except pandas.computation.ops.UndefinedVariableError:
			die("Error: annotation file filter specified variable that does not exist, filter was: %s" % annot_filter,logger)

		if df_anno is None or df_anno.shape[0] == 0:
			die("Error: annotation file filter '%s' resulted in an empty data set" % annot_filter,logger)

	return df_anno

def create_group_file(outfile,df_anno,gene_col,logger=None):
	"""
	Given an annotation file, create a group file.
	:param outfile: Output filepath for group file
	:param df_anno: Annotation data frame (should have already been filtered if requested by user)
	:param gene_col: Name of the gene column
	:param logger: Optional logger object
	:return:
	"""

	if logger is not None:
		logger.debug("called create_group_file: \n" + pprint.pformat(locals()))

	# Write out group file
	with open(outfile,"w") as out:
		for gene, group in df_anno.groupby(gene_col):
			variants = group["EPACTS"].unique()

			print >> out, "\t".join([gene] + list(variants))

def check_file_exists(filename):
	if os.path.exists(filename):
		return True
	else:
		print >> sys.stderr, "Error: file " + filename + " does not exist\n"
		sys.exit(1)

def die(msg,logger=None):
	if logger is not None:
		logger.error(msg)
	else:
		print >> sys.stderr, msg

	sys.exit(1)