#!/usr/bin/env python
DESCRIPTION = '''
Predict homomerization of ColabFold-AlphaFold2 structure

NOTE:
 - Can be given either the output directory from colabfold_batch (will pick top ranked 
   relaxed or unrelaxed structure) or a specific PDB file to analyze.
 - The PDB file used needs to have been generated from a di-oligomeric prediction.
   That is, ColabFold needs to have been run using two copies of the protein (ABC:ABC)
   So that symmetry can be estimated.
'''
import sys
import os
import shutil
import argparse
import logging
import gzip
from pathlib import Path



## Pass arguments.
def main():
	## Pass command line arguments. 
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
	parser.add_argument('results',
		type=str, 
		help='One of: 1) directory with colabfold_batch PDB files, or 2) a specific PDB file'
	)
	parser.add_argument('oligomerization', 
		type=str, 
		help='Oligomerization results output directory.'
	)
	parser.add_argument('--debug', 
		required=False, action='store_true', 
		help='Print DEBUG info (default: %(default)s)'
	)
	args = parser.parse_args()
	
	## Set up basic debugger
	logFormat = "[%(levelname)s]: %(message)s"
	logging.basicConfig(format=logFormat, stream=sys.stderr, level=logging.INFO)
	if args.debug:
		logging.getLogger().setLevel(logging.DEBUG)
	
	logging.debug('%s', args) ## DEBUG
	
	colabfold_predict_homomerization(args.results, args.oligomerization)



def colabfold_predict_homomerization(results, oligomerization):
	# Create output directory
	os.makedirs(oligomerization, exist_ok=True)
	
	# For each PDB file that we have found
	for prefix, pdb_file in get_pdb_file(results):
		prefix = os.path.basename(prefix) # Remove directory name from prefix
		logging.info(f"Processing PDB file: {pdb_file}; Output Prefix: {prefix}") ## INFO
		
		# Output variables
		ini_pdb = f"{oligomerization}/{prefix}.pdb"
		sym_pdb = f"{oligomerization}/{prefix}_symm.pdb"
		sym_def = f"{oligomerization}/{prefix}_symm.txt"
		sym_log = f"{oligomerization}/{prefix}_symm.log.txt"
		sym_mon = f"{oligomerization}/{prefix}_symm.monomer_count.txt"
		
		# Copy delected PDB file to output dir
		shutil.copyfile(pdb_file, ini_pdb)
		
		# Run 'make_symmdef_file.pl' to predict symmetry
		#  - FROM: https://files.ipd.uw.edu/krypton/make_symmdef_file.pl
		make_symmdef_file = Path(os.path.dirname(os.path.realpath(__file__)) + '/make_symmdef_file.pl')
		if not make_symmdef_file.exists():
			raise OSError(f"Cant find make_symmdef_file.pl, please download from https://files.ipd.uw.edu/krypton/make_symmdef_file.pl")
		cmd = f"perl {make_symmdef_file} -m NCS -p {ini_pdb} 2> {sym_log} 1> {sym_def}"
		logging.debug(f"Running cmd: {cmd}") ## DEBUG
		ret = os.system(cmd)
		
		# Print make_symmdef_file.pl run log
		for line in open(sym_log, "r").readlines():
			print(line.rstrip())
		
		# Fail if command returned non-zero error code.
		if ret != 0:
			raise Exception(f'make_symmdef_file.pl failed to finish correctly! See the above log for details.')
		
		# Get the number of copies
		copies = 0
		for line in open(sym_pdb, "r").readlines():
			if line.startswith('TER'):
				copies += 1
		
		# Write monomer count to file
		#logging.info(f"Estimated a total of {copies} monomers in the symmetric complex.") ## INFO
		with open(sym_mon,"w") as f:
			f.write(f"{copies}\n")
		
		# Print explanation of output files
		logging.info(f"Best selected structure PDB file: {ini_pdb}") ## INFO
		logging.info(f"Infered symmetric structure PDB file: {sym_pdb}") ## INFO
		logging.info(f"Number of estimated monomers in comples: {sym_mon}") ## INFO
		logging.info(f"Log file: {sym_log}") ## INFO
		logging.info(f"Defintion file: {sym_def}") ## INFO



def get_pdb_file(input_path):
	'''
	Read one of: 1) directory with colabfold_batch PDB files, or 2) a specific PDB file
	'''
	
	# Check 'input_path' exists in some form
	input_path = Path(input_path)
	if not input_path.exists():
		raise OSError(f"{input_path} could not be found")
	
	## Get PDB file from 'input_path'
	# 'input_path' is a file
	if input_path.is_file():
		if input_path.suffix == ".pdb":
			pdb_file = input_path
			prefix   = str(input_path).rstrip(".pdb")
			yield(prefix, pdb_file)
		else:
			raise ValueError(f"Unknown file format {input_path.suffix}")
	
	# 'input_path' is a dir
	else:
		assert input_path.is_dir(), "Expected either an input file or a input directory"
		
		# For each A3M file in 'input_path' (infer the number of queries)
		prefixes = []
		for file in sorted(input_path.iterdir()):
			if not file.is_file():
				continue
			if file.suffix.lower() not in [".a3m"]:
				continue
			prefixes.append(str(file).rstrip(".a3m"))
		if len(prefixes) == 0:
			raise ValueError(f"No *.a3m files found in {input_path}")
		
		logging.debug(f"A3M prefixes found: {prefixes}") ## DEBUG
		
		# For each A3M prefix (i.e., each protein query in results)
		for prefix in prefixes:
			queries = []
			for file in sorted(input_path.iterdir()):
				logging.debug(f"Checking {prefix} {file}") ## DEBUG
				if not file.is_file():
					continue
				if file.suffix.lower() not in [".pdb"]:
					continue
				if not str(file).startswith(prefix):
					continue
				queries.append(file)
			
			# Get either top ranked relaxed structure, or if relaxation not performed, top ranked unrelazed structure
			for file in queries:
				if '_relaxed_rank_001_' in str(file):
					pdb_file = file
					break
				if '_unrelaxed_rank_001_' in str(file):
					pdb_file = file
					break
			else:
				raise ValueError(f"No top relaxed or unrelaxed PDB file found in {queries}")
			
			# Yield each set
			yield(prefix, pdb_file)



if __name__ == '__main__':
	main()
