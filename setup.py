from setuptools import setup

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
      ],
      packages=['MutationInfo'],
)


