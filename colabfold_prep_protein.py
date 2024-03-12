#!/usr/bin/env python
DESCRIPTION = '''
Prepare protein sequence for input into ColabFold-AlphaFold2

Will cleanup protein sequence of non-standard characters
 - remove terminal * from protein
 - replace X and internal * with 'G'
    - Non-standard amino acids will crash alphafold so need 
      to replace them with a neutral non-polar amino acid

Will also produce proteins formatted for homo-oligomer prediction if required.
 - This script can only prepare homo-oligomers, protein files for hetero-oligomers will
   need to be prepared manually.

E.G.,
 - If '-n=4', the following sequences
>1
ABCD
>2
ABC
 - Will be convered to
>1
ABCD:ABCD:ABCD:ABCD
>2
ABC:ABC:ABC:ABC
 - This will result in ColabFold predicting each sequence as a quad homo-oligomer
'''
import sys
import os
import argparse
import logging
import gzip
from itertools import groupby
import re

## Pass arguments.
def main():
	## Pass command line arguments. 
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
	parser.add_argument('infasta',
		type=lambda x: File(x, 'r'), 
		help='Input [gzip] file'
	)
	parser.add_argument('outfasta', 
		type=lambda x: File(x, 'w'), 
		help='Output [gzip] file'
	)
	parser.add_argument('-n', '--num_oligomers', 
		required=False, default=1, type=int,
		help='Format proteins so ColabFold predictes a homo-oligomer of size N (default: %(default)s)'
	)
	parser.add_argument('-q', '--quiet',
		required=False, action='store_true',
		help='Dont print info about each sequence with characters removed (default: %(default)s)'
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
	elif args.quiet:
		logging.getLogger().setLevel(logging.WARNING)
	
	logging.debug('%s', args) ## DEBUG
	
	
	with args.infasta as infasta, args.outfasta as outfasta:
		colabfold_prep_protein(infasta, outfasta, args.num_oligomers)


def colabfold_prep_protein(infasta, outfasta, num_oligomers, ambiguous_char = "G", allowed_chars = set(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])):
	for header, seq in fasta_iter(infasta):
		# Remove stop codon
		seq = re.sub('\\*$', '', seq)
		
		# Identify non-standard characters and remove any if found
		seq_set = set(seq.upper())
		char_diff = seq_set.difference(allowed_chars)
		if len(char_diff) > 0:
			logging.info('Non-standard characters %s found in sequence %s', char_diff, header) ## DEBUG
			for c in char_diff:
				seq = seq.replace(c, ambiguous_char)
		
		# Duplicate cleaned sequence for homo-oligomer prediction
		seq = ':'.join([seq] * num_oligomers)
		
		# Write seq
		outfasta.write(f">{header}\n{seq}\n")


class File(object):
	'''
	Context Manager class for opening stdin/stdout/normal/gzip files.

	 - Will check that file exists if mode='r'
	 - Will open using either normal open() or gzip.open() if *.gz extension detected.
	 - Designed to be handled by a 'with' statement (other wise __enter__() method wont 
	    be run and the file handle wont be returned)
	
	NOTE:
		- Can't use .close() directly on this class unless you uncomment the close() method
		- Can't use this class with a 'for' loop unless you uncomment the __iter__() method
			- In this case you should also uncomment the close() method as a 'for'
			   loop does not automatically cloase files, so you will have to do this 
			   manually.
		- __iter__() and close() are commented out by default as it is better to use a 'with' 
		   statement instead as it will automatically close files when finished/an exception 
		   occures. 
		- Without __iter__() and close() this object will return an error when directly closed 
		   or you attempt to use it with a 'for' loop. This is to force the use of a 'with' 
		   statement instead. 
	
	Code based off of context manager tutorial from: https://book.pythontips.com/en/latest/context_managers.html
	'''
	def __init__(self, file_name, mode):
		## Upon initializing class open file (using gzip if needed)
		self.file_name = file_name
		self.mode = mode
		
		## Check file exists if mode='r'
		if not os.path.exists(self.file_name) and mode == 'r':
			raise argparse.ArgumentTypeError("The file %s does not exist!" % self.file_name)
	
		## Open with gzip if it has the *.gz extension, else open normally (including stdin)
		try:
			if self.file_name.endswith(".gz"):
				#print "Opening gzip compressed file (mode: %s): %s" % (self.mode, self.file_name) ## DEBUG
				self.file_obj = gzip.open(self.file_name, self.mode+'b')
			else:
				#print "Opening normal file (mode: %s): %s" % (self.mode, self.file_name) ## DEBUG
				self.file_obj = open(self.file_name, self.mode)
		except IOError as e:
			raise argparse.ArgumentTypeError('%s' % e)
	def __enter__(self):
		## Run When 'with' statement uses this class.
		#print "__enter__: %s" % (self.file_name) ## DEBUG
		return self.file_obj
	def __exit__(self, type, value, traceback):
		## Run when 'with' statement is done with object. Either because file has been exhausted, we are done writing, or an error has been encountered.
		#print "__exit__: %s" % (self.file_name) ## DEBUG
		self.file_obj.close()
#	def __iter__(self):
#		## iter method need for class to work with 'for' loops
#		#print "__iter__: %s" % (self.file_name) ## DEBUG
#		return self.file_obj
#	def close(self):
#		## method to call .close() directly on object.
#		#print "close: %s" % (self.file_name) ## DEBUG
#		self.file_obj.close()



def fasta_iter(fh):
    	"""
    	Given a fasta file. yield tuples of header, sequence
	Clears description from seq name.
	
	From: https://www.biostars.org/p/710/
	Updated: 11/09/2018
	Version: 0.2
	"""
    	# ditch the boolean (x[0]) and just keep the header or sequence since
    	# we know they alternate.
    	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    	for header in faiter:
        	# drop the ">" and description
        	header = next(header)[1:].strip().split(' ')[0]
        	# join all sequence lines to one.
        	seq = "".join(s.strip() for s in next(faiter))
        	yield header, seq


if __name__ == '__main__':
	main()
