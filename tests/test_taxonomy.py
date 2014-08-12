#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "2014, The Viroscape Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"
 
from numpy import array
from cogent.util.unit_test import TestCase, main
from biom.table import DenseTable
from taxonomy import biom_table_from_blast_results,split_taxonomy_dict_on_semicolons

class TaxonomyTests(TestCase):
    """ """
    
    def setUp(self):
        self.rap_result_with_taxonomy =\
          [RAP_HEADER_WITH_TAXONOMY]
        
        self.rap_result_with_taxonomy.extend(RAP_DATA_LINES_WITH_TAXONOMY)
        
        self.exp_biom_table_from_rap_with_taxonomy =\
          EXP_BIOM_TABLE_RAP_LINES_WITH_TAXONOMY         
    
    def split_taxonomy_dict_on_semicolons_valid_input(self):
        """split_taxonomy_dict_on_semicolons functions with valid input"""
        input_data ={\
          "obs_id1":{'taxonomy': 'Eukaryota    ;Metazoa;Chordata  ;Mammalia;Chiroptera;Vespertilionidae;Scotophilus;Scotophilus kuhlii','other_metadata':"don't;split;this"},\
          "obs_id2":{'taxonomy':"Eukaryota;Metazoa;Chordata;Aves;Anseriformes;Anatidae;Cairina;Cairina moschata",'other_metadata':"don't;split;this;either"}}

        obs = split_taxonomy_dict_on_semicolons(input_data)
        exp = {"obs_id1":{'taxonomy':['Eukaryota','Metazoa','Chordata','Mammalia','Chiroptera','Vespertilionidae','Scotophilus','Scotophilus kuhlii'],'other_metadata':"don't;split;this"},"obs_id2":{'taxonomy': ['Eukaryota','Metazoa','Chordata','Aves','Anseriformes','Anatidae','Cairina','Cairina moschata'],'other_metadata':"don't;split;this;either"}}
        self.assertEqualItems(obs,exp)

EXP_BIOM_TABLE_RAP_LINES_WITH_TAXONOMY =\
"""# Constructed from biom file\n#OTU ID\tlane7-index02-CGATGT-N3.txt_filtered\tlane7-index03-TTAGGC-N4.txt_filtered\ttaxonomy
11866\t1.0\t1.0\t['Viruses', 'Unknown', 'Unknown', 'Unknown', 'Unknown', 'Retroviridae', 'Alpharetrovirus', 'Avian_myeloblastosis_virus']
32630\t1.0\t0.0\t['Unknown', 'Unknown', 'Unknown', 'Unknown', 'Unknown', 'Unknown', 'Unknown', 'synthetic___construct']
9598\t1.0\t0.0\t['Eukaryota', 'Metazoa', 'Chordata', 'Mammalia', 'Primates', 'Hominidae', 'Pan', 'Pan_troglodytes']
9606\t0.0\t1.0\t['Eukaryota', 'Metazoa', 'Chordata', 'Mammalia', 'Primates', 'Hominidae', 'Homo', 'Homo_sapiens']"""

RAP_HEADER_WITH_TAXONOMY = "\t".join(["# Fields:","Query","Subject","identity","aln-len","mismatch","gap-openings","q.start","q.end","s.start","s.end","log(e-value)","bit-score","sampleid","gi","taxid","taxonomy"])

RAP_DATA_LINES_WITH_TAXONOMY = [\
"\t".join(["DB775P1:233:D1RBTACXX:7:1103:10381:4391","gi|307686499|dbj|BAJ21180.1|","100","32","0","0","2","97","261","292","-11.73","76.6406","lane7-index02-CGATGT-N3.txt_filtered","307686499","32630","Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;synthetic   construct"]),\
"\t".join(["DB775P1:233:D1RBTACXX:7:1107:13823:85667","gi|114625270|ref|XP_520079.2|","100","32","0","0","3","98","1622","1653","-12.66","79.7221","lane7-index02-CGATGT-N3.txt_filtered","114625270","9598","Eukaryota;Metazoa;Chordata;Mammalia;Primates;Hominidae;Pan;Pan troglodytes"]),\
"\t".join(['DB775P1:233:D1RBTACXX:7:1106:5991:43041', 'gi|127590|sp|P01104.2|MYB_AVIMB', '93.9394', '33', '2', '0', '1', '99', '42', '74', '-11.38', '75.485', 'lane7-index02-CGATGT-N3.txt_filtered', '127590', '11866', 'Viruses;Unknown;Unknown;Unknown;Unknown;Retroviridae;Alpharetrovirus;Avian myeloblastosis virus']),\
"\t".join(['DB775P1:233:D1RBTACXX:7:1101:15166:39365', 'gi|765157|gb|AAB31928.1|', '100', '32', '0', '0', '3', '98', '624', '655', '-11.26', '75.0998', 'lane7-index03-TTAGGC-N4.txt_filtered', '765157', '11866', 'Viruses;Unknown;Unknown;Unknown;Unknown;Retroviridae;Alpharetrovirus;Avian myeloblastosis virus']),\
"\t".join(['DB775P1:233:D1RBTACXX:7:1109:18474:85187', 'gi|2398621|emb|CAA04019.1|', '100', '33', '0', '0', '1', '99', '1360', '1392', '-11.26', '75.0998', 'lane7-index03-TTAGGC-N4.txt_filtered', '2398621', '9606', 'Eukaryota;Metazoa;Chordata;Mammalia;Primates;Hominidae;Homo;Homo sapiens'])
]


if __name__ == "__main__":
    main()
