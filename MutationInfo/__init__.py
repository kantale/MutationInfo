
import os
import json

from Bio import Entrez
from appdirs import *

#Entrez.email = global_vars['email']

__docformat__ = 'reStructuredText'

#  http://thomas-cokelaer.info/tutorials/sphinx/docstring_python.html 

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
			if not _directory_exists(local_directory):
				raise EnvironmentError('Local Directory %s does not exist' % (str(local_directory)))
			self.local_directory = local_directory

		self._properties_file = os.path.join(self.local_directory, self._properties_file)
		if not self._file_exists(self._properties_file):
			#Create property file
			with open(self._properties_file, 'w') as f:
				f.write('{}\n')

		#Read property file
		self.properties = self._load_json_filename(self._properties_file)

		#Get email
		if not email is None:
			self.properties['email'] = email
		elif not 'email' in self.properties:
				self.properties['email'] = raw_input('I need an email to query Entrez. Please insert one: ')
		Entrez.email = self.properties['email']
		print 'Using email for accessing Entrez: %s' % (str(Entrez.email))

		#Save properties file
		self._save_json_filenane(self._properties_file, self.properties)

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
		if not self._file_exists(filename):
			with open(filename, 'w') as f:
				f.write(fasta)

	def _load_ncbi_fasta_filename(self, ncbi_access_id):
		'''
		Load NCBI fasta file
		'''
		filename = self._ncbi_fasta_filename(self, ncbi_access_id)
		if not self._file_exists(filename):
			return None

		with open(filename) as f:
			fasta = f.read()

		return fasta

	@staticmethod
	def _directory_exists(dirname):
		'''
		Check if directory exists
		'''
		return os.path.isdir(dirname)

	@staticmethod
	def _file_exists(filename):
		'''
		Check if filename exists
		'''
		return os.path.isfile(filename) 

	@staticmethod
	def _mkdir_p(dirname):
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
	def _get_application_dir():
		'''
		Create a cross platform local directory for this app
		Reference: https://pypi.python.org/pypi/appdirs/1.4.0 
		'''
		directory = user_data_dir('MutationInfo', '')
		if not MutationInfo._directory_exists(directory):
			MutationInfo._mkdir_p(directory)

		return directory

	@staticmethod
	def _load_json_filename(filename):
		'''
		Load a json file
		'''
		with open(filename) as f:
			data = json.load(f)

		return data

	@staticmethod
	def _save_json_filenane(filename, data):
		with open(filename, 'w') as f:
			f.write(json.dumps(data, indent=4) + '\n')

