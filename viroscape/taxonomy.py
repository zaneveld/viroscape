#!/usr/bin/env python
# File created on 2 Jul 2014
from __future__ import division
from collections import defaultdict

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

def split_taxonomy_dict_on_semicolons(taxonomy_dict,target_md_field="taxonomy"):
    """Return a list of stripped entries produced by splitting taxa_string on semicolons"""
    result = defaultdict(dict)
    #Keys are typically observation_ids
    #Values are dicts of metadata taxonomy strings

    #Split each value on semicolons and strip 
    #to generate a dict of lists
    for k,v in taxonomy_dict.iteritems():
        fields = v[target_md_field].split(';')
        result[k][target_md_field]=[f.strip().replace(" ","_") for f in fields]
    return dict(result)

def yield_blast_results_with_taxonomy(results,header_map,ncbi_taxonomy_tree,gi_to_taxid_mapping):
    """Yield successive lines of RAP or BLAST results , annotated with taxonomy.
    
    results -- m8 format similarity search results (only RAP tested currently)
    header_map -- dict mapping column names in results to their index
    ncbi_taxonomy_tree -- PhyloNode object representation of the NCBI
      taxonomy tree
    gi_to_taxid_mapping -- a mapping between NCBI gis and taxids 
    
    
    """
    subject_idx = header_map["Subject"]
    assigned_gis = 0
    unassigned_gis = 0
    for fields in results:
        last_field = len(fields)
        subject_id = fields[subject_idx]  #
        gi = parse_gi_from_subject_id(subject_id)
        taxid = gi_to_taxid_mapping.get(gi,"Unknown")
        if taxid == "Unknown":
            unassigned_gis +=1
            continue
        else:
            assigned_gis +=1
        lineage = get_lineage_from_taxid(int(taxid.strip()),ncbi_taxonomy_tree)
        fields.extend([f.strip() for f in [gi,taxid,lineage]])
        yield fields
    print "total unassigned gis: %i" %unassigned_gis
    print "total assigned gis: %i" %assigned_gis


def gis_from_blast_lines(lines):
    """Extract subject gi numbers from blast results"""
    header_line,header_map,blast_result =\
      parse_similarity_search_lines(lines)
    print header_line,header_map,blast_result
    all_gis = set()
    for fields in blast_result:
        subject = fields[header_map["Subject"]]
        all_gis.add(parse_gi_from_subject_id(subject))
    print "gis_from_blast_lines: number of gis:",len(all_gis)
    return all_gis



def get_lineage_from_taxid(taxid,ncbi_taxonomy_tree,ranks=\
  ['superkingdom','kingdom','phylum','class','order','family','genus','species']):
    """Return a lineage for a given taxid
    taxid -- an NCBI taxon id as an integer
    ncbi_taxonomy_tree -- a tree object for the NCBI taxonomy 
      (from PyCogent's NcbiTaxonomyFromFiles)
    
    
    """
    try:
        node = ncbi_taxonomy_tree.ById[taxid]
    except KeyError:
        #Couldn't find taxid in tree
        print "Couldn't find taxid '%s' in taxonomy tree. Assigning 'Unknown'." %(taxid)
        return ";".join(["Unknown"]*len(ranks))
    lineage = get_lineage(node, ranks)
    return lineag



def get_lineage(node, ranks):
    """Return a semicolon-delimited taxonomic lineage given a node and list of ranks

    node -- a node object in the NCBI taxonomic tree
    ranks -- a list of ranks to output. For example ["family","species"] 

    """
    #Make a dict of rank name by index in the list of ranks
    ranks_lookup = dict([(r,idx) for idx,r in enumerate(ranks)])

    #Start by setting all ranks to Unknown
    lineage = ["Unknown"] * len(ranks)

    curr = node
    while curr.Parent is not None:
        if curr.Rank in ranks_lookup:
            lineage[ranks_lookup[curr.Rank]] = curr.Name
        curr = curr.Parent

    return ";".join(map(str,lineage))

