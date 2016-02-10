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

* If the variant is in rs format, then 
    * Access the [UCSC tables](https://genome.ucsc.edu/cgi-bin/hgTables) through the [cruzdb](https://pypi.python.org/pypi/cruzdb) package. 
    * If this fails then try the [Variant Effect Predictor](http://asia.ensembl.org/Tools/VEP) through the [pyVEP](https://github.com/kantale/pyVEP) package.
* If the variant is in HGVS then:
    * Try to parse the variant with the [biocommon/hgvs](https://bitbucket.org/biocommons/hgvs) parser. 
    * If the parse fails then look if the variant contains some common mistakes in HGVS formatting. Correct if possible and then try again. For example remove parenthesis in the following variant: `NM_001042351.1:-1923(A>C)`
    * If parse still fails then make a request to the [mutalyzer.nl](https://mutalyzer.nl/). For example `NT_005120.15:c.IVS1-72T>G` is parsed only from mutalyzer but not from biocommons/hgvs
    * If biocommons/hgvs parses the variant then use the [variantmapper](http://hgvs.readthedocs.org/en/latest/examples/manuscript-example.html#project-genomic-variant-to-a-new-transcript) method to locate the location of the variant in the reference assembly.
    * If this method fails (for example `M61857.1:c.121A>G` crashes mutalyzer!) then use the [pyhgvs](https://github.com/counsyl/hgvs) package and the `hgvs_to_vcf` method to convert the variant in a [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) entry.
    * If this method fails then look at the [LOVD](http://www.lovd.nl/3.0/home) database.

If all the aforementioned methods fail then: 

* Download the FASTA sequence of the trascript of the variant from [NCBI database](http://www.ncbi.nlm.nih.gov/nuccore).
* If the position of the variant is in coding (c.) coordinates then convert to genomic (g.) coordinates. To do that, we use the [Coordinate mapper](https://github.com/lennax/biopython/tree/f_loc5/Bio/SeqUtils/Mapper) addition of biopython.
* Perform a [blat search](https://genome.ucsc.edu/cgi-bin/hgBlat?command=start) from UCSC. This methods performs an alignment search of the fasta sequence in the reference assembly. In case this succeeds then report the location of the variant in the reference genome. 

## Installation 
To install MutationInfo, download the package and run:
```bash
python setup.py install
```
Important! Requires 13 GB of disk space.

## Contact 
[Alexandros Kanterakis](mailto:kantale@ics.forth.gr)

