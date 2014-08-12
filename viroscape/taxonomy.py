#!/usr/bin/env python
# File created on 2 Jul 2014
from __future__ import division


__author__ = "Jesse Zaneveld"
__copyright__ = "Vega Thurber Lab Copyright 2014"
__credits__ = ["Jesse Zaneveld",""] 
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

"""
Library functions for viral taxonomy
"""
def add_taxonomy_to_blast_results(lines,ncbi_taxonomy_tree,gi_to_taxid_mapping):
    """Output an annotated BLAST result with fields for taxid & lineage
    
    lines -- blast of a similarity search result 
       NOTE: the subject id is assumed to contain the gi in a format
        like this: gi|148686800|gb|EDL18747.1|
    
    header_map -- a dictionary mapping the names of columns to their 
      index for example if column 0 is 'Query', header_map['Query']=0.
    
    ncbi_taxonomy_tree -- a PyCogent tree object representing 
     the NCBI taxonomy tree made with PyCogents NCBITaxonomyTreeFromFiles,
      and NCBI's nodes.dmp and names.dmp files.
    
    gi_to_taxid_mapping -- a mapping between the relevant gis and taxids 
     
    """
    #The Subject field will always be in the same column
    #Look up ahead of time which column that is
    header_line,header_map,blast_result =\
      parse_similarity_search_lines(lines)

    #Make a generator object that will yield parsed BLAST results 
    results_iterator =\
      yield_blast_results_with_taxonomy(blast_result,header_map,\
      ncbi_taxonomy_tree,gi_to_taxid_mapping)

    #Update the header map and header line with new entries for 'gi',
    #'taxid',and 'taxonomy'
    
    last_field_idx = max(header_map.values())
    header_map["gi"] = last_field_idx +1
    header_map["taxid"] = last_field_idx +2
    header_map["taxonomy"] = last_field_idx +3
    new_header_line = header_line.strip()
    print "old header line:",header_line
    new_header_line = "\t".join([new_header_line,"gi","taxid","taxonomy"])
    print "new header line:",new_header_line
    return new_header_line,header_map,results_iterator
