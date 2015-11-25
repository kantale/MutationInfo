
import os
import re
import json
import logging

from Bio import Entrez
from appdirs import *

#Entrez.email = global_vars['email']

__docformat__ = 'reStructuredText'

"""
TODO: 
* More documentation   http://thomas-cokelaer.info/tutorials/sphinx/docstring_python.html 
* Fix setup.py http://stackoverflow.com/questions/3472430/how-can-i-make-setuptools-install-a-package-thats-not-on-pypi 
"""

class MutationInfo(object):
	"""The MutationInfo class contains methods to get variant information 

	"""


	_properties_file = 'properties.json'

	def __init__(self, local_directory=None, email=None):

		#Get local directory
		if local_directory is None:
			self.local_directory = self._get_application_dir()
			print 'Using local directory: %s' % (self.local_directory) 
		else:
			if not Utils.directory_exists(local_directory):
				raise EnvironmentError('Local Directory %s does not exist' % (str(local_directory)))
			self.local_directory = local_directory

		self._properties_file = os.path.join(self.local_directory, self._properties_file)
		if not Utils.file_exists(self._properties_file):
			#Create property file
			with open(self._properties_file, 'w') as f:
				f.write('{}\n')

		#Read property file
		self.properties = Utils.load_json_filename(self._properties_file)

		#Get email
		if not email is None:
			self.properties['email'] = email
		elif not 'email' in self.properties:
				self.properties['email'] = raw_input('I need an email to query Entrez. Please insert one: ')
		Entrez.email = self.properties['email']
		print 'Using email for accessing Entrez: %s' % (str(Entrez.email))

		#Save properties file
		Utils.save_json_filenane(self._properties_file, self.properties)

	@staticmethod
	def fuzzy_hgvs_corrector(variant, transcript=None, ref_type=None):
		"""
		Try to correct a wrong HGVS-ish variant by checking if it matches some patterns with common mistakes.
		Following directions from here: http://www.hgvs.org/mutnomen/recs-DNA.html#sub
		This is by far not exhaustive.. 

		:param variant: The name of the variant (example: 1234A>G)
		:param trascript: In case the variant does not have a transcript part then use this.
		:param ref_type: In case the variant does not include a reference type indicator (c or g) the define it here
		"""

		if ref_type not in [None, 'c', 'g']:
			raise ValueError('Available values for ref_type: None, "c" and "g" . Found: %s' % (str(ref_type)))

		#Exclude variants in unicode
		new_variant = str(variant)

		#Check if we have all necessary information
		if not ':' in new_variant:
			if transcript is None:
				logging.error('Variant: %s does not include a transcript part (":") and the transcript argument is None. Returning None ' % (new_variant))
				return None

			search = re.search(r'[cg]\.', new_variant)
			if search is None:
				if ref_type is None:
					logging.error('Variant: %s does not include a reference type part (c or g) and the ref_type argument is None. Returning None ' % (new_variant))
					return None
			new_variant = str(transcript) + ':' + ref_type + '.' + new_variant

		#Case 1
		#Instead if ">" the input is: "->". For example: 
		if '->' in new_variant:
			logging.warning('Variant: %s  . "->" found. Substituting it with ">"' % (new_variant))
			new_variant = new_variant.replace('->', '>')

		# Case 2
		# The variant contains / in order to declare two possible substitutions 
		search = re.search(r'([ACGT])>([ACGT])/([ACGT])', new_variant)
		if search:
			logging.warning('Variant: %s  . "/" found suggesting that this contains 2 variants' % (new_variant))
			new_variant_1 = re.sub(r'([ACGT])>([ACGT])/([ACGT])', r'\1>\2', new_variant)
			new_variant_2 = re.sub(r'([ACGT])>([ACGT])/([ACGT])', r'\1>\3', new_variant)
			return [fuzzy_hgvs_corrector(new_variant_1), fuzzy_hgvs_corrector(new_variant_2)]

		return new_variant


	def get_info(variant):
		"""

		:param variant: A NCBI access ID (example: 'NM_000463.2')

		"""


		
		pass

	def _get_fasta_from_nucleotide_entrez(self, ncbi_access_id, pure=False): # For example NG_000004.3
		'''
		handle4 = Entrez.efetch(db='nucleotide', id='101011606', rettype='fasta', retmode='text') 

		'''
		
		fasta = self._load_ncbi_fasta_filename(ncbi_access_id)
		if fasta is None:

			print 'Could not find local file for %s Querying Entrez..' % (ncbi_access_id)
			handle = Entrez.efetch(db='nuccore', id=ncbi_access_id, retmode='text', rettype='fasta')
			fasta = handle.read()

			if pure:
				#Remove initial tags and newlines
				fasta = ''.join(fasta.split('\n')[1:])

			print 'Fasta fetched: %s...' % (fasta[0:20])
		
			self._save_ncbi_fasta_filename(ncbi_access_id, fasta)

		return fasta


	def _ncbi_fasta_filename(self, ncbi_access_id):
		'''
		Create filename that contains NCBI fasta file
		'''
		return os.path.join(self.local_directory, ncbi_access_id + '.fasta')

	def _save_ncbi_fasta_filename(self, ncbi_access_id, fasta):
		'''
		Save NCBI fasta to file
		'''

		filename = self._ncbi_fasta_filename(self, ncbi_access_id)
		if not Utils.file_exists(filename):
			with open(filename, 'w') as f:
				f.write(fasta)

	def _load_ncbi_fasta_filename(self, ncbi_access_id):
		'''
		Load NCBI fasta file
		'''
		filename = self._ncbi_fasta_filename(self, ncbi_access_id)
		if not Utils.file_exists(filename):
			return None

		with open(filename) as f:
			fasta = f.read()

		return fasta

	@staticmethod
	def _get_application_dir():
		'''
		Create a cross platform local directory for this app
		Reference: https://pypi.python.org/pypi/appdirs/1.4.0 
		'''
		directory = user_data_dir('MutationInfo', '')
		if not Utils.directory_exists(directory):
			Utils.mkdir_p(directory)

		return directory



class Utils(object):
	'''
	Useful functions to help manange files
	'''

	@staticmethod
	def directory_exists(dirname):
		'''
		Check if directory exists
		'''
		return os.path.isdir(dirname)



	@staticmethod
	def file_exists(filename):
		'''
		Check if filename exists
		'''
		return os.path.isfile(filename) 


	@staticmethod
	def mkdir_p(dirname):
		'''
		Create directory
		Functionality similar with: mkdir -p
		Reference: http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python 
		'''
		try:
			os.makedirs(dirname)
		except OSError as exc: # Python >2.5
			if exc.errno == errno.EEXIST and os.path.isdir(dirname):
				pass
			else:
				raise

	@staticmethod
	def load_json_filename(filename):
		'''
		Load a json file
		'''
		with open(filename) as f:
			data = json.load(f)

		return data

	@staticmethod
	def save_json_filenane(filename, data):
		'''
		Save a json file
		'''
		with open(filename, 'w') as f:
			f.write(json.dumps(data, indent=4) + '\n')

def test():
	'''
	Testing cases
	'''

	print MutationInfo.fuzzy_hgvs_corrector('1048G->C')
	print MutationInfo.fuzzy_hgvs_corrector('1048G->C', transcript='NM_001042351.1')
	try: 
		MutationInfo.fuzzy_hgvs_corrector('1048G->C', transcript='NM_001042351.1', ref_type='p')
	except Exception as e:
		print 'Exception:', str(e)
	print MutationInfo.fuzzy_hgvs_corrector('1048G->C', transcript='NM_001042351.1', ref_type='c')
