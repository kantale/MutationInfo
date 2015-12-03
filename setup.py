
import sys
import subprocess

try:
	from setuptools import setup
except ImportError as e:
	print 'setuptools package is not installed'
	print 'On linux install with  wget https://bootstrap.pypa.io/ez_setup.py -O - | python '
	print 'For more info please visit: https://pypi.python.org/pypi/setuptools'
	sys.exit(1)


try:
	import Bio
except ImportError as e:
	print 'Biopython is not installed'
	# Check ig gcc exists.. 
	try:
		subprocess.call(['gcc', '--version'])
	except OSError as e:
		print '***** WARNING ******'
		if e.errno == os.errno.ENOENT:
			print 'Before installing biopython you may need to install a C compiler such as gcc'
			print 'On Linux you can install with the following command:'
			print 'sudo apt-get install gcc python-dev libpq-dev python-pip '
			
		else:
			print 'Trying to run gcc produced the following error: %s' % (str(e))
			print 'Not being able to run gcc might be a problem when installing biopython'


setup(name='MutationInfo',
      version='0.0.1',
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
            'pyhgvs>=2.0.0beta',
      ],
      dependency_links=['https://github.com/counsyl/hgvs/tarball/master#egg=pyhgvs-2.0.0beta',],
      packages=['MutationInfo'],
)


