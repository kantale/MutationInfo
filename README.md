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



Although there are numerous tools that extract meta-information of variants, these tools are highly specific. This tools mainly is a wrapper to 


Requires 13 GB disk space.

