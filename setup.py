
import os
import re
import sys
import subprocess

#Check python version
if sys.version_info[:2] >= (3,0):
      print 'Sorry.. MutationInfo is a python 2 tool..'
      sys.exit(1)

if sys.version_info[:2] < (2.7):
      print 'Python 2.7 is required for MutationInfo'
      print 'Python version detected: ', str(sys.version_info)
      sys.exit(1)

try:
	from setuptools import setup
except ImportError as e:
	print 'setuptools package is not installed'
	print 'On linux install with the following command:'
	print 'wget https://bootstrap.pypa.io/ez_setup.py -O - | sudo python '
	print 'For more info please visit: https://pypi.python.org/pypi/setuptools'
	sys.exit(1)

try:
	import Bio
except ImportError as e:
      print 'Biopython is not installed'
      # Check if gcc exists.. 
      try:
            subprocess.call(['gcc', '--version'])
      except OSError as e:
            print '***** WARNING ******'
            if e.errno == os.errno.ENOENT:
                  print 'Before installing biopython you may need to install a C compiler such as gcc.'
                  print 'gcc and other necessary libraries can be installed on Linux with the following command:'
                  print 'sudo apt-get install gcc python-dev libpq-dev python-pip python-mysqldb-dbg'			
            else:
                  print 'Trying to run gcc produced the following error: %s' % (str(e))
                  print 'Not being able to run gcc might be a problem when installing biopython'
            print '********************'
            _ = raw_input("Press Enter to continue or Ctrl-C to stop.. ")


try:
      #subprocess.call(['pg_config', '--version'])
      pg_v = None
      pg_version = subprocess.check_output(['pg_config', '--version']) # i.e. 'PostgreSQL 9.4.4\n' 
      pg_s = re.search(r'([\d]+)\.([\d]+)\.([\d]+)', pg_version)
      pg_v = (int(pg_s.group(1)), int(pg_s.group(2)), int(pg_s.group(3)))

except:
      print 'PostgreSQL does not seem to be installed'
      print 'PostgreSQL is a dependency for BioPython (http://biopython.org/DIST/docs/biosql/python_biosql_basic.html#htoc2)'
      _ = raw_input("Press Enter to continue or Ctrl-C to stop.. ")

# psycopg2 does not seem to work for before a certain version of PostgreSQL 
# We took the last working version from here: https://github.com/psycopg/psycopg2/blob/master/tests/test_connection.py#L434
# Also: https://github.com/kantale/MutationInfo/issues/16 
if pg_v:
      if pg_v < (9,2,0):
            print 'Your current PostgreSQL version is: %s' % (str(pg_v))
            print 'This is lower than 9.2.0 which is the highest version of PostgreSQL that is compatible with psycopg2.'
            print 'psycopg2 is the PostgreSQL driver that certain MutationInfo tools (biocommons hgvs) are using.'
            _ = raw_input("Press Enter to continue or Ctrl-C to stop.. ")



setup(name='MutationInfo',
      version='1.1.0',
      description='Tool to retrieve meta-information of genetic variants',
      url='https://github.com/kantale/MutationInfo',
      author='Alexandros Kanterakis',
      author_email='alexandros.kanterakis@gmail.com',
      license='MIT',
      # https://pypi.python.org/pypi?%3Aaction=list_classifiers 
      classifiers=[
            'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      install_requires=[
            'biopython',
            'appdirs',
            'hgvs>=0.4,<0.5',
            'feedparser',
            'cruzdb',
            'pygr',
            'sqlalchemy',
            'beautifulsoup4',
            'pyhgvs>=0.0.1',
            'pyVEP>=0.0.1',
            'myvariant', # https://pypi.python.org/pypi/myvariant/0.2.0 
      ],
      # http://stackoverflow.com/questions/3472430/how-can-i-make-setuptools-install-a-package-thats-not-on-pypi 
      dependency_links=[
            'https://github.com/counsyl/hgvs/tarball/master#egg=pyhgvs-2.0.0',
            'https://github.com/kantale/pyVEP/tarball/master#egg=pyVEP-2.0.0',
      ],
      packages=['MutationInfo', 'biopython_mapper'],
)

# Check if psycopg2 is 'importable'
try:
      import psycopg2
except ImportError as e:
      if 'Library not loaded: libssl.1.0.0.dylib' in str(e):
            print '='*10 + '==========' + '='*10
            print ' '*10 + 'IMPORTANT:'
            print '='*10 + '==========' + '='*10
            print 'Module psycopg2 although installed cannot be imported properly. Error message:'
            print '=' * 20
            print str(e)
            print '=' * 20
            print 'To resolve this, before running MutationInfo set the following environment variable:'
            lib_path = os.path.split(os.path.split(os.__file__)[0])[0]
            DYLD_FALLBACK_LIBRARY_PATH = os.environ.get('DYLD_FALLBACK_LIBRARY_PATH', '') # Not used..
            command = "export DYLD_FALLBACK_LIBRARY_PATH={}:$DYLD_FALLBACK_LIBRARY_PATH".format(lib_path)
            print command
            print 'For more please check: http://stackoverflow.com/questions/27264574/import-psycopg2-library-not-loaded-libssl-1-0-0-dylib'

      raise e


try:
      import hgvs
except ImportError as e:
      if 'cannot import name ExtendedInterpolation' in str(e):
            print str(e)
            print 'This is a known issue.'
            print 'Please refer to https://github.com/kantale/MutationInfo/issues/9 in order to resolve it'
      raise e

