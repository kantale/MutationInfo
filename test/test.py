
import unittest

import logging
logging.basicConfig(level=logging.DEBUG)

from MutationInfo import MutationInfo

mi = MutationInfo()

class TestMutationInfo(unittest.TestCase):

    maxDiff = None

    def test_FYZZY_HGVS_CORRECTOR(self):
        print '------FUZZY HGVS CORRECTOR---------'
        ret = MutationInfo.fuzzy_hgvs_corrector('1048G->C')
        print ret
        self.assertEqual(ret, None)

        ret = MutationInfo.fuzzy_hgvs_corrector('1048G->C', transcript='NM_001042351.1')
        print ret
        self.assertEqual(ret, None)

        try: 
            MutationInfo.fuzzy_hgvs_corrector('1048G->C', transcript='NM_001042351.1', ref_type='p')
        except Exception as e:
            self.assertEqual(str(e), 'Available values for ref_type: None, "c" and "g" . Found: p')

        ret = MutationInfo.fuzzy_hgvs_corrector('1048G->C', transcript='NM_001042351.1', ref_type='c')
        print ret
        self.assertEqual(ret, 'NM_001042351.1:c.1048G>C')

        ret = MutationInfo.fuzzy_hgvs_corrector('1387C->T/A', transcript='NM_001042351.1', ref_type='c')
        print ret
        self.assertEqual(ret, ['NM_001042351.1:c.1387C>T', 'NM_001042351.1:c.1387C>A'])

        ret = MutationInfo.fuzzy_hgvs_corrector('1387C->T/A')
        print ret
        self.assertEqual(ret, None)

        ret = MutationInfo.fuzzy_hgvs_corrector('-1923(A>C)', transcript='NT_005120.15', ref_type='g')
        print ret
        self.assertEqual(ret, 'NT_005120.15:g.-1923A>C')

        ret = MutationInfo.fuzzy_hgvs_corrector('NT_005120.15:c.1160(CC>GT)')
        print ret
        self.assertEqual(ret, 'NT_005120.15:c.1160_1161delinsGT')


    def test_HGVS_PARSER(self):
        print '--------HGVS PARSER-----------------'
        ret = MutationInfo.biocommons_parse('unparsable')
        print ret
        self.assertEqual(ret, None)

    def test_MUTALYZER(self):
        print '-------Mutalyzer---------------------'
        ret = mi._search_mutalyzer('NT_005120.15:c.IVS1-72T>G', gene='UGT1A1')
        print ret
        self.assertEqual(ret, 'NT_005120.15:g.608362T>G')


    def test_LOVD(self):
        print '--------LOVD------------------------'
        ret = mi.lovd_transcript_dict['NM_000367.2']
        print ret
        self.assertEqual(ret, [u'TPMT', u'hg19'])

        ret = mi._search_lovd('NM_000367.2', 'c.-178C>T')
        print ret
        self.assertEqual(ret, (u'6', 18155397, None, u'hg19'))

    def test_VARIATION_REPORTER(self):
        print '--------VARIATION REPORTER------------------------'
        ret = mi.get_info('NM_099999999.2:c.166C>T', method='VARIATION_REPORTER')
        print ret
        self.assertIsNone(ret)

        ret = mi.get_info('NM_099999.2:c.166C>T', method='VARIATION_REPORTER')
        print ret
        self.assertIsNone(ret)

        ret = mi.get_info('NM_017781.2:c.166C>T', method='VARIATION_REPORTER')
        print ret
        self.assertEqual(ret, {'chrom': '7', 'notes': 'Variation Reporter did c. to g conversion', 'source': 'NC_transcript', 'genome': 'GRCh37.p13', 'offset': 1023013, 'alt': 'T', 'ref': 'C'})

    def test_TRANSVAR(self):
        print '--------TRANSVAR------------------------'
        ret = mi.get_info('NM_017781.2:c.166C>T', method='TRANSVAR')
        print ret
        self.assertEqual(ret, {'chrom': '7', 'notes': 'Transvar conerted NM_017781.2:c.166C>T to chr7:g.1023013C>T', 'source': 'counsyl_hgvs_to_vcf', 'genome': 'hg19', 'offset': 1023013, 'alt': 'T', 'ref': 'C'})

    def test_GET_INFO_HGVS(self):
        print '--------GET INFO HGVS--------------------'

        info = mi.get_info('1387C->T/A', transcript='NM_001042351.1', ref_type='c') # This should try first to correct 
        print info
        self.assertEqual(info,  [{'chrom': 'x', 'notes': '', 'source': 'NC_transcript', 'genome': 'GRCh37.p13', 'offset': 153760473, 'alt': 'A', 'ref': 'G'}, {'chrom': 'x', 'notes': '', 'source': 'NC_transcript', 'genome': 'GRCh37.p13', 'offset': 153760473, 'alt': 'T', 'ref': 'G'}])

        info = mi.get_info('1387C->T/A') # This should fail (return None)
        print info
        self.assertEqual(info,  None)

        info =  mi.get_info('NM_006446.4:c.1198T>G')
        print info
        self.assertEqual(info,  {'chrom': '12', 'notes': '', 'source': 'NC_transcript', 'genome': 'GRCh37.p13', 'offset': 21355487, 'alt': 'G', 'ref': 'T'})
        
        info = mi.get_info('XYZ_006446.4:c.1198T>G')
        print info
        self.assertEqual(info,  None)
        
        info = mi.get_info('NM_006446.4:c.456345635T>G')
        print info
        self.assertEqual(info,  {'chrom': '12', 'notes': 'Variant: NM_006446.4:c.456345635T>G . biocommons method splign reported error: The given coordinate is outside the bounds of the reference sequence. , '
    'Variant: NM_006446.4:c.456345635T>G . biocommons method blat reported error: The given coordinate is outside the bounds of the reference sequence. , '
    'Variant: NM_006446.4:c.456345635T>G . biocommons method genewise method failed: No alignments for NM_006446.4 in GRCh37 using genewise', 'source': 'Mutalyzer', 'genome': 'GRCh38.p7', 'offset': 477582748, 'alt': 'G', 'ref': 'T'})

        info = mi.get_info('NG_000004.3:g.253133T>C')
        print info
        self.assertEqual(info, {'chrom': '7', 'notes': u'MUTALYZER POSITION CONVERTER REPORTED ERROR: The Accession number NG_000004 could not be found in our database (or is not a chromosome). / MUTALYZER POSITION CONVERTER REPORTED ERROR: The Accession number NG_000004 could not be found in our database (or is not a chromosome). / Mutalyzer did c_g conversion, INFO BY BLAT', 'source': 'Mutalyzer', 'genome': 'hg19', 'offset': 99264573, 'alt': 'G', 'ref': 'A'})

        info = mi.get_info({})
        print info
        self.assertEqual(info, None)

        info = mi.get_info(['NM_001042351.1:c.1387C>T', 'NM_001042351.1:c.1387C>A'])
        print info
        self.assertEqual(info, [{'chrom': 'x', 'notes': '', 'source': 'NC_transcript', 'genome': 'GRCh37.p13', 'offset': 153760473, 'alt': 'A', 'ref': 'G'}, {'chrom': 'x', 'notes': '', 'source': 'NC_transcript', 'genome': 'GRCh37.p13', 'offset': 153760473, 'alt': 'T', 'ref': 'G'}])

        info = mi.get_info('NC_000001.11:g.97593343C>A')
        print info
        self.assertEqual(info, {'chrom': '1', 'notes': '', 'source': 'NC_transcript', 'genome': 'GRCh38.p7', 'offset': 97593343, 'alt': 'A', 'ref': 'C'})

        info = mi.get_info('M61857.1:c.121A>G') # No exons in genbank file
        print info
        info_to_check = {k:v for k,v in info.iteritems() if k != 'notes'}
        self.assertEqual(info_to_check,  {'chrom': '10', 'source': 'BLAT', 'genome': 'hg19', 'offset': 96698560, 'alt': 'G', 'ref': 'A'})

        info = mi.get_info('J02843.1:c.-1295G>C') # Test biopython c2g
        print info
        self.assertEqual(info,  {'chrom': '10', 'notes': u'Variant: J02843.1:c.-1295G>C . biocommons method splign method failed: No alignments for J02843.1 in GRCh37 using splign / Variant: J02843.1:c.-1295G>C . biocommons method blat method failed: No alignments for J02843.1 in GRCh37 using blat / Variant: J02843.1:c.-1295G>C . biocommons method genewise method failed: No alignments for J02843.1 in GRCh37 using genewise / MUTALYZER POSITION CONVERTER REPORTED ERROR: The accession number J02843 could not be found in our database (or is not a transcript). / MUTALYZER POSITION CONVERTER REPORTED ERROR: The accession number J02843 could not be found in our database (or is not a transcript). / Variant: J02843.1:c.-1295G>C . Mutalyzer returned the following critical error: No gene specified. Please choose from: ', 'source': 'BLAT', 'genome': 'hg19', 'offset': 135339602, 'alt': 'C', 'ref': 'G'})

        info = mi.get_info('J02843.1:c.7632T>A') # Test biopython c2g, this position is out of exon boundaries
        print info
        self.assertEqual(info,  None)
        
        info = mi.get_info('X17059.1:c.1091insAAA') # No exons and no mRNA. This is meant to fail
        print info
        self.assertEqual(info,  None)

        info = mi.get_info('AY545216.1:g.8326_8334dupGTGCCCACT')
        print info
        self.assertEqual(info,  {'chrom': '22', 'notes': u'MUTALYZER POSITION CONVERTER REPORTED ERROR: The Accession number AY545216 could not be found in our database (or is not a chromosome). / MUTALYZER POSITION CONVERTER REPORTED ERROR: The Accession number AY545216 could not be found in our database (or is not a chromosome). / Mutalyzer did c_g conversion, INFO BY BLAT', 'source': 'Mutalyzer', 'genome': 'hg19', 'offset': 42522668, 'alt': '', 'ref': 'CACGGGTGA'})

        info = mi.get_info('NT_005120.15:c.-1126C>T', gene='UGT1A1')
        print info
        self.assertEqual(info,  {'chrom': '2', 'notes': u'Variant: NT_005120.15:c.-1126C>T . biocommons method splign method failed: No alignments for NT_005120.15 in GRCh37 using splign / Variant: NT_005120.15:c.-1126C>T . biocommons method blat method failed: No alignments for NT_005120.15 in GRCh37 using blat / Variant: NT_005120.15:c.-1126C>T . biocommons method genewise method failed: No alignments for NT_005120.15 in GRCh37 using genewise / MUTALYZER POSITION CONVERTER REPORTED ERROR: The accession number NT_005120 could not be found in our database (or is not a transcript). / MUTALYZER POSITION CONVERTER REPORTED ERROR: The accession number NT_005120 could not be found in our database (or is not a transcript). / Variant: NT_005120.15(UGT1A1):c.-1126C>T . Mutalyzer returned the following critical error: C not found at position 600562, found G instead. / Variant: NT_005120.15:c.-1126C>T . ***SERIOUS*** Reference on fasta (A) and Reference on variant name (C) are different!', 'source': 'BLAT', 'genome': 'hg19', 'offset': 234069238, 'alt': 'T', 'ref': 'C'})

        info = mi.get_info('NM_000367.2:c.-178C>T')
        print info
        self.assertEqual(info,  {'chrom': '6', 'notes': '', 'source': 'counsyl_hgvs_to_vcf', 'genome': 'hg19', 'offset': 18155397, 'alt': 'A', 'ref': 'G'})

        info = mi.get_info('NT_005120.15:c.IVS1-72T>G', gene='UGT1A1')
        print info
        self.assertEqual(info,  {'chrom': '2', 'notes': u'MUTALYZER POSITION CONVERTER REPORTED ERROR: The Accession number NT_005120 could not be found in our database (or is not a chromosome). / MUTALYZER POSITION CONVERTER REPORTED ERROR: The Accession number NT_005120 could not be found in our database (or is not a chromosome). / Mutalyzer did c_g conversion, INFO BY BLAT', 'source': 'Mutalyzer', 'genome': 'hg19', 'offset': 234675608, 'alt': 'G', 'ref': 'T'})


    def test_GET_INFO_RS(self):
        print '--------GET INFO RS--------------------'
        # Testing rs SNPs
        info = mi.get_info('rs53576')
        print info
        self.assertEqual(info,  {'chrom': '3', 'notes': '', 'source': 'VEP', 'genome': u'GRCh38', 'offset': 8762685, 'alt': u'G', 'ref': u'A'})

        info = mi.get_info('rs4646438') # insertion variation 
        print info
        self.assertEqual(info,  {'chrom': '7', 'notes': '', 'source': 'VEP', 'genome': u'GRCh38', 'offset': 99766411, 'alt': u'T', 'ref': u''})

        info = mi.get_info('rs305974') # This SNP is not in UCSC 
        print info
        self.assertEqual(info,  {'chrom': '19', 'notes': '', 'source': 'VEP', 'genome': u'GRCh38', 'offset': 41121963, 'alt': u'A', 'ref': u'G'})

        info = mi.get_info('rs773790593') # Both UCSC and VEP fail
        print info
        self.assertEqual(info,  {'chrom': '22', 'notes': '', 'source': 'VEP', 'genome': u'GRCh38', 'offset': 42130778, 'alt': u'A', 'ref': u'G'})

        info = mi.get_info('rs758320086') # Deletion detected by VEP
        print info
        self.assertEqual(info,  {'chrom': '22', 'source': 'VEP', 'genome': u'GRCh38', 'offset': 42128249, 'alt': u'', 'ref': u'AGTT',  'notes': ''})

        info = mi.get_info('rs1799752') # Deletion from UCSC , # UCSC Returnes lengthTooLong 
        print info
        self.assertEqual(info,  {'chrom': '17', 'notes': '', 'source': 'VEP', 'genome': u'GRCh38', 'offset': 63488529, 'alt': [u'ATACAGTCACTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCC', u'G'], 'ref': u''})

        info = mi.get_info('rs72466463') # Insertion from UCSC
        print info
        self.assertEqual(info,  {'chrom': '2', 'notes': '', 'source': 'VEP', 'genome': u'GRCh38', 'offset': 38071144, 'alt': u'GGTGGCATGA', 'ref': u''})

        info = mi.get_info('rs8175347') # UCSC short tandem repeat (microsatellite) variation 
        print info
        self.assertEqual(info,  {'chrom': '2', 'notes': '', 'source': 'VEP', 'genome': u'GRCh38', 'offset': 233760236, 'alt': [u'TATATATATATATATA', u'TATATATATATATA', u'TATATATATATA'], 'ref': [u'(TA)5', u'(TA)6', u'(TA)7', u'(TA)8']})

        info = mi.get_info('rs28399445') # UCSC Substitution
        print info
        self.assertEqual(info,  {'chrom': '19', 'notes': '', 'source': 'VEP', 'genome': u'GRCh38', 'offset': 40848265, 'alt': u'T', 'ref': u'GC'}) #  Check https://mutalyzer.nl/name-checker?description=NM_000762.5%3Ac.608_609delGCinsA 

        info = mi.get_info('rs267607275') # UCSC Multiple alleles 
        print info
        self.assertEqual(info,  {'chrom': '6', 'notes': '', 'source': 'VEP', 'genome': u'GRCh38', 'offset': 18149126, 'alt': [u'C', u'G'], 'ref': u'A'} )

        info = mi.get_info('NG_008377.1:g.6502_6507delCTCTCT') # BLAT Deletion 
        print info
        self.assertEqual(info,  {'chrom': '19', 'notes': u'MUTALYZER POSITION CONVERTER REPORTED ERROR: The Accession number NG_008377 could not be found in our database (or is not a chromosome). / MUTALYZER POSITION CONVERTER REPORTED ERROR: The Accession number NG_008377 could not be found in our database (or is not a chromosome). / Mutalyzer did c_g conversion, INFO BY BLAT', 'source': 'Mutalyzer', 'genome': 'hg19', 'offset': 41354851, 'alt': '', 'ref': 'GAGAGA'})

        info = mi.get_info('rs113940699') # Fails VEP, succeeds on MyVariantInfo
        print info
        self.assertEqual(info, {'chrom': '22', 'notes': 'Variant: rs113940699 . VEP returned an empty list', 'source': 'MyVariantInfo', 'genome': 'hg19', 'offset': 42522748, 'alt': u'T', 'ref': u'C'})

        info = mi.get_info('rs72549356') # Fails all three! 
        print info
        self.assertIsNone(info)

if __name__ == '__main__':
    '''
    Run: 
        * python test.py , to run all tests
        * python test.py TestMutationInfo.test_VARIATION_REPORTER , to run a specific test: http://stackoverflow.com/questions/15971735/running-single-test-from-unittest-testcase-via-command-line
        * python test.py TestMutationInfo.test_TRANSVAR
    '''

    unittest.main()
