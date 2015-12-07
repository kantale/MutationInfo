# MutationInfo

MutationInfo is a python package to extract the position, the reference and the alternative sequence of a genomic variant. It accepts variants in [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/) rs format or in [HGVS](http://www.hgvs.org/mutnomen/) format. 

The main purpose of MutationInfo is to simplify the process of locating a variant in a dataset (i.e. of sequences or variants) that is aligned in a human reference genome (for example hg19 or hg28). It mainly wraps a collection of existing tools with a simple interface.

Example:
```python
from MutationInfo import MutationInfo
mi = MutationInfo()

mi.get_info('rs53576')
{'ref': 'A', 'alt': 'G', 'chrom': 'chr3', 'genome': 'hg19', 'offset': 8804371L}

mi.get_info('NM_006446.4:c.1198T>G')
{'ref': 'T', 'alt': 'G', 'chrom': '12', 'genome': 'GRCh37.p13', 'offset': 21355487}
```

## How it works 

MutationInfo tries to infer the position, reference and alternative of a variant through the following pipeline:

* If the variant is in rs format, then access the [UCSC tables](https://genome.ucsc.edu/cgi-bin/hgTables) through the [cruzdb](https://pypi.python.org/pypi/cruzdb) package. 
* If the variant is in HGVS then:
    * Try to parse the variant with the [biocommon/hgvs](https://bitbucket.org/biocommons/hgvs) parser. 
    * If the parse fails then look if the variant contains some common mistakes in HGVS formatting. Correct if possible and then try again. For example remove parenthesis in the following variant: `NM_001042351.1:-1923(A>C)`
    * If parse still fails then make a request to the [mutalyzer.nl](https://mutalyzer.nl/). For example `NT_005120.15:c.IVS1-72T>G` is parsed only from mutalyzer and not from biocommons/hgvs
    

there are numerous tools that extract meta-information of variants, these tools are highly specific. This tools mainly is a wrapper to 


Requires 13 GB disk space.

