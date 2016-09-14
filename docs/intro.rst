Introduction
============

The main purpose of MutationInfo is to simplify the process of locating a variant in a dataset (i.e. of sequences or variants) 
that is aligned in a human reference genome (for example hg19 or hg38).
It mainly wraps a collection of existing tools with a simple interface.

Example:

.. code-block:: python

  from MutationInfo import MutationInfo
  mi = MutationInfo()

  # RS variant
  print mi.get_info('rs53576')
  {'chrom': '3', 'notes': '', 'source': 'UCSC',
  'genome': 'hg19', 'offset': 8804371L, 'alt': 'G', 'ref': 'A'}

  # HGVS variant
  print mi.get_info('NM_000367.2:c.-178C>T')
  {'chrom': '6', 'notes': '', 'source': 'counsyl_hgvs_to_vcf', 
  'genome': 'hg19', 'offset': 18155397, 'alt': 'A', 'ref': 'G'}

How it works 
------------

MutationInfo tries to infer the position, reference and alternative of a variant through the following pipeline:

* If the variant is in rs format, then
    * Try the `Variant Effect Predictor <http://asia.ensembl.org/Tools/VEP>`_ through the `pyVEP <https://github.com/kantale/pyVEP>`_ package.
    * If this fails, try the `MyVariant.info <http://myvariant.info/>`_ service.
    * If this fails, access the `UCSC tables <https://genome.ucsc.edu/cgi-bin/hgTables>_ through the `cruzdb <https://pypi.python.org/pypi/cruzdb>`_ package. 
* If the variant is in HGVS then:
    * Try to parse the variant with the `biocommon/hgvs <https://bitbucket.org/biocommons/hgvs>`_ parser. 
    * If the parse fails then look if the variant contains some common mistakes in HGVS formatting. Correct if possible and then try again. For example remove parenthesis in the following variant: ``NM_001042351.1:-1923(A>C)``
    * If parse still fails then make a request to the `mutalyzer.nl <https://mutalyzer.nl/>`_ . For example ``NT_005120.15:c.IVS1-72T>G`` is parsed only from mutalyzer but not from biocommons/hgvs
    * If neither of these methods are able to parse the variant then return ``None``. 
    * If biocommons/hgvs parses the variant then use the `variantmapper <http://hgvs.readthedocs.org/en/latest/examples/manuscript-example.html#project-genomic-variant-to-a-new-transcript>`_ method to locate the location of the variant in the reference assembly.
    * If this method fails then use the `pyhgvs <https://github.com/counsyl/hgvs>`_ package and the `hgvs_to_vcf` method to convert the variant in a `VCF <https://en.wikipedia.org/wiki/Variant_Call_Format>`_ entry.
    * If this method fails then use Mutalyzer's `Position Converter <https://mutalyzer.nl/position-converter>`_ 
    * if this method fails then use the Mutalyzer's `Name Checker <https://mutalyzer.nl/>`_ which generates a genomic description of the variant. Then perform a `blat search <https://genome.ucsc.edu/cgi-bin/hgBlat?command=start>`_) with this variant (see below).
    * If both methods from Mutalyzer fail (for example ``M61857.1:c.121A>G`` crashes mutalyzer!) then:
        * Download the FASTA sequence of the trascript of the variant from `NCBI database <http://www.ncbi.nlm.nih.gov/nuccore>`_.
        * If the position of the variant is in coding (c.) coordinates then convert to genomic (g.) coordinates. To do that, we use the `Coordinate mapper <https://github.com/lennax/biopython/tree/f_loc5/Bio/SeqUtils/Mapper>`_ addition of biopython.
        * Perform a `blat search <https://genome.ucsc.edu/cgi-bin/hgBlat?command=start>`_) from UCSC. This methods performs an alignment search of the fasta sequence in the reference assembly. In case this succeeds then report the location of the variant in the reference genome. 
    * If this method fails then search the `LOVD <http://www.lovd.nl/3.0/home>`_ database.
    * If all the aforementioned methods fail then return ``None``



Installation 
============

.. note::
  Important! Requires 13 GB of disk space.

To install MutationInfo, download the latest release from https://github.com/kantale/MutationInfo/releases , uncompress and run:

.. code-block:: bash

  python setup.py install


Then the first time you instantiate the MutationInfo class, it installs all required datasets:

.. code-block:: python

  from MutationInfo import MutationInfo
  mi = MutationInfo()

Installation in Ubuntu
----------------------

Before installing in Ubuntu Linux, make sure that the following packages / tools are installed:

.. code-block:: bash

  sudo apt-get update
  sudo apt-get install git
  sudo apt-get install gcc python-dev libpq-dev python-pip python-mysqldb-dbg

  wget https://bootstrap.pypa.io/ez_setup.py -O - | sudo python


Test Installation
-----------------

To verify that everything works fine run: ``python test.py`` in ``test/`` directory. The output after the long log messages should be:

::

  ----------------------------------------------------------------------
  Ran 6 tests in 21.923s

  OK

Troubleshooting
---------------

Possible problems from installing / running MutationInfo are:

* Exception: ``psycopg2.OperationalError: invalid connection option "application_name"``
   See also: https://github.com/kantale/MutationInfo/issues/16 . Most likely, the version of PostgreSQL in your system is too old. 
* Exception: ``ImportError: cannot import name ExtendedInterpolation``
   See also: https://github.com/kantale/MutationInfo/issues/9 . One solution is to downgrade the ``future`` package. In that case, it is a good practice to 
   run MutationInfo in a virtualenv so that the whole system is not affected.
* Exception: ``ImportError: No module named MySQLdb`` 
   See also: https://github.com/kantale/MutationInfo/issues/7 . mysql is not installed in the system.
* Error Message: ``Library not loaded: libssl.1.0.0.dylib`` 
   See: https://github.com/kantale/MutationInfo/issues/5 . 


