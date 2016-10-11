[![DOI](https://zenodo.org/badge/45920824.svg)](https://zenodo.org/badge/latestdoi/45920824)

# MutationInfo

MutationInfo is a python package to extract the position, the reference and the alternative sequence of a genomic variant. It accepts variants in [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/) rs format or in [HGVS](http://www.hgvs.org/mutnomen/) format. 

The main purpose of MutationInfo is to simplify the process of locating a variant in a dataset (i.e. of sequences or variants) that is aligned in a human reference genome (for example hg19 or hg38). It mainly wraps a collection of existing tools with a simple interface.

Example:
```python
from MutationInfo import MutationInfo
mi = MutationInfo()

mi.get_info('rs53576')
{'chrom': '3', 'source': 'UCSC', 'genome': 'hg19', 
 'offset': 8804371L, 'alt': 'G', 'ref': 'A'}

mi.get_info('NM_006446.4:c.1198T>G')
{'chrom': '12', 'source': 'NC_transcript', 'genome': 'GRCh37.p13', 
 'offset': 21355487, 'alt': 'G', 'ref': 'T'}
```

# Documentation 
The documentation is here: http://mutationinfo.readthedocs.io/en/latest/ 

## License 
MIT License (MIT)

## Contact 
[Alexandros Kanterakis](mailto:kantale@ics.forth.gr)

