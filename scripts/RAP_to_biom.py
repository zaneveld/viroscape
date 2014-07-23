#!/usr/bin/env python
# File created on 2 Jul 2014
from __future__ import division


#of credits
__author__ = "Jesse Zaneveld"
__copyright__ = "Vega Thurber Lab Copyright 2014"
__credits__ = ["Jesse Zaneveld","Stephanie Rosales","RyanMcMinds"] 
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = ""
__email__ = ""
__status__ = "Development"
 

from cogent.util.option_parsing import parse_command_line_parameters, make_option
from cogent.parse.blast import MinimalBlastParser9
from cogent.parse.ncbi_taxonomy import NcbiTaxonomyFromFiles
from cogent.util.misc import create_dir
from collections import defaultdict
from os import path,listdir
from os.path import isdir,splitext,exists
from biom.table import table_factory,SparseOTUTable
from numpy import zeros

from viroscape.parse import parse_delimited_data_lines,parse_similarity_search_lines,parse_gi_from_subject_id,make_header_map

script_info = {}
script_info['brief_description'] = "This script reads a list of gi numbers (each can appear more than once), and generates a BIOM table based on NCBI taxonomy"
script_info['script_description'] = "For each query sequence"
script_info['script_usage'] = [("","< examples of script usage go here>. Predict metagenome from count file and OTU table.","%prog -i blast_results_m8.txt --e_value 1e-5 --percent_similarity 70  -o filtered_blast_results_m8.txt")]
script_info['output_description']= "Outputs a table of similarities filtered by user specified values."
script_info['required_options'] = [
 make_option('-i','--input_blast_results_file',help='Input BLAST or RAP search results file, in tabular (-m 8) format. Alternatively, a directory can be supplied. In that case, all files in the directory must be -m 8 format RAP search results'),
 make_option('-o','--output',type="new_filepath",help='output directory for results')
]
script_info['optional_options']=[
 make_option('-e','--e_value',action="store",type="float",default=1e-10,help='The maximal e-value of the records that will be preserved. [default: %default]'),
 make_option('-p','--percent_aligned',type="float",default=0.70,help='The maximal e-value of the records that will be preserved. [default: %default]'),
 make_option('--gi_to_taxid_file',default="/raid2/labs/Thurber_lab/databases/ncbi_taxonomy/gi_taxid_prot.dmp ",
  help='A file mapping gi numbers to NCBI taxon ids (be careful to use the protein or nucleotide .dmp file as appropriate). [default: %default]'),
 make_option('--ncbi_taxonomy_nodes_file',default="/raid2/labs/Thurber_lab/databases/ncbi_taxonomy/nodes.dmp",\
  help="Path to the NCBI taxonomy nodes.dmp file [default:%default]"),\
 make_option('--ncbi_taxonomy_names_file',default="/raid2/labs/Thurber_lab/databases/ncbi_taxonomy/names.dmp",\
  help="Path to the NCBI taxonomy names.dmp file [default:%default]"),\
 make_option('--filename_to_sampleid_map_fp',default=None,\
  help="Path to a tab-delimited file mapping the names of the input blast files to sampleids [default:%default]"),\
 make_option('--recover',action='store_true',default=False,\
  help="Assume a previous run failed.  Recover one step back from the failure  [default:%default]")

]
script_info['version'] = __version__

def biom_table_from_blast_results(lines,
  sample_column_name='sampleid',observation_column_name='taxid',\
  metadata_cols=['taxonomy'],delimiter="\t"):
    """Build a BIOM table object from a tab-delimited file with observation, 
    sample, and 1+ metadata columns

    lines -- tab-delimited BLAST or RAP search result lines
    observation_column_name -- the header label for the observation id col
    sample_column_name -- the label for the sampleid column
    metadata_cols -- a list of labels for columns that will be used as metadata
    delimiter -- the delimiter separating columns
    """
    #Need to have BLAST lines already annotated with sample id and taxonomy
    #Step 1. 1st pass create a set of sample ids and observation ids 
    
    sample_ids,observation_ids,observation_metadata =\
      biom_data_from_blast_results(lines,sample_column_name,\
      observation_column_name)
    
    unique_sample_ids = sorted(list(set(sample_ids)))
    unique_observation_ids = sorted(list(set(observation_ids)))
    
    sample_map = {sample_id:i for i,sample_id in enumerate(unique_sample_ids)}
    observation_map = {observation_id:i for i,observation_id in enumerate(unique_observation_ids)}
    data_as_dict = defaultdict(int)
      
    for sample_id,observation_id in zip(sample_ids,observation_ids):
        data_as_dict[(observation_map[observation_id],sample_map[sample_id])] += 1
    
    data_array = array_from_dict(data_as_dict)
    print "Obs length:",len(observation_metadata)
    print "Data.shape:",data_array.shape
    print "Contruting biom table..."
    result = table_factory(data_array,unique_sample_ids,\
      unique_observation_ids,\
      constructor=SparseOTUTable)

    print "Adding observation metadata to biom table"
    #Oddly, the table_factory constructor in biom doesn't accept dict
    #input, so observation metadata is added as a separate method call
    
    observation_metadata_fields =\
      split_taxonomy_dict_on_semicolons(observation_metadata)
    
    result.addObservationMetadata(observation_metadata_fields)
    return result

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

def array_from_dict(sparse_dict):
    """Convert a dict of {(row,col):value} to a dense array"""
    #The size of the array should be one larger
    #than the largest index
    N = max(r for r,c in sparse_dict.iterkeys())+1
    M = max(c for r,c in sparse_dict.iterkeys())+1

    result = zeros((N,M))
    for coords,counts in sparse_dict.iteritems():
        row,col = coords
        result[row][col]=counts
    return result

def biom_data_from_blast_results(lines,\
  sample_column_name='sampleid',observation_column_name='taxid',metadata_cols=['taxonomy']):    
    """Retrieve sample ids, observations and metadata from blast or RAP result lines""" 
    
    
    sample_ids = []
    observation_ids = []
    metadata={} 
    header_line,header_map,data_fields =\
      parse_similarity_search_lines(lines,strict=True) 
    sample_idx = header_map[sample_column_name]
    observation_idx = header_map[observation_column_name]
    metadata_indices = [(m,header_map[m]) for m in metadata_cols]
    print "Metadata indices:",metadata_indices
    for data_line_number,fields in enumerate(data_fields):
        fields = [f for f in fields]
        #print "fields:",fields
        #print "Data line #%i" % data_line_number
        #print "i,Fields:",[(i,curr_f) for i,curr_f in enumerate(fields)]
        #print "fields:",fields
        #print "len(fields):",len(fields)
        #print "header_map:",header_map
        #print "sample_idx:",sample_idx
        #print "observation_idx:",observation_idx
        observation_id = fields[observation_idx]
        sample_id = fields[sample_idx]
        observation_ids.append(observation_id)
        sample_ids.append(sample_id)
        
        #Build a dict of metadata for this row only
        curr_metadata = {}
        for metadata_label,i in metadata_indices:
            metadata_field = fields[i]
            curr_metadata[metadata_label]=metadata_field
        
        #If an observation id occurs more than once,
        #we simply overwrite the metadata.
        #All observation metadata *should* be identical for 
        #each instance of the same observation id.
        #It is presently the responsibility of the user
        #to ensure this is the case (no check is performed)
        metadata[observation_id]=curr_metadata
    print "Metadata returned by biom_table_from_blast_lines:",metadata  
    return sample_ids,observation_ids,metadata
 
        

def best_RAP_results(lines,limit_to_sampleid=False,large_e_value=10000.0):
    """Return a dictionary of field indices, and a dict of the best RAP results by query"""
    
    header_line,header_map,blast_results=parse_similarity_search_lines(lines)
    
    #If filtering by sampleid, make sure its present
    if 'sampleid' not in header_map.keys():
        raise ValueError("best_RAP_results cannot filter by sampleid because 'sampleid' column is not present in header_map: %s" %header_map)
    
    query_field_idx = header_map['Query']
    log_e_val_idx = header_map['log(e-value)']
    #Define a function that will return a large float by default
    #when the defaultdict (below) is queried for a missing value
    def large_float_factory():
        return large_e_value

    e_vals_by_query = defaultdict(large_float_factory)
    result_by_query = {}
    
    for i,blast_result in enumerate(blast_results):
        if not blast_result:
            continue
        if limit_to_sampleid:
            sampleid = blast_result[header_map['sampleid']]
            #Skip records not from the sample of interest
            if sampleid != limit_to_sampleid:
                continue
        #Get some key values form the list of results by index
        
        query_id = blast_result[query_field_idx]
        log_e_value = blast_result[log_e_val_idx]
        
        #Convert log e-values to e-values
        e_value = log_e_value_to_e_value(log_e_value)

        if e_value >=  e_vals_by_query[query_id]:
            #Not a better BLAST hit than whats there
            continue

        #Update the current best e_val and id
        #to be this id
        e_vals_by_query[query_id] = e_value
        result_by_query[query_id] = blast_result
    
    
    return header_line,header_map,result_by_query.values()                     

def log_e_value_to_e_value(log_e_value,base=10.0):
    l = float(log_e_value)
    e_value = base**l
    return e_value 

def filter_RAP_results(lines,max_e_val,min_percent_aligned, gi_only=False,\
    add_sampleid=False, include_header=True):    
    """Filter RAP search (m8 format) results, yielding results"""
    #blast_results = MinimalBlastParser9(input_file)
    header_line,header_map,blast_results =\
      parse_similarity_search_lines(lines)
    
    
    if add_sampleid:
        #Add a column to the header for the sampleid
        header_line = "\t".join([header_line.strip(),'sampleid'])+"\n"
        # Assign this new column a key in the header map
        last_field_idx = max(header_map.values())
        sampleid_col_idx = last_field_idx + 1
        header_map['sampleid']=sampleid_col_idx

    if include_header:
        yield [header_line]

    #Before parsing the full results, find which field indexes are 
    #the ones we care about
    log_e_val_index = header_map['log(e-value)']
    percent_aligned_index = header_map['identity'] 
    subject_index = header_map['Subject']
    #These files can be large (many GB) so we never want to construct
    #a list of all results in memory
    for i,blast_result in enumerate(blast_results):
        try:
            log_e_value = float(blast_result[log_e_val_index])
        except IndexError,e:
            print "BAD BLAST RESULT:",blast_result
            print "log_e_val_index:",log_e_val_index
            raise IndexError(e)
        
        percent_aligned = float(blast_result[percent_aligned_index]) 
        e_value = log_e_value_to_e_value(log_e_value)
        
        if e_value > float(max_e_val):
            continue 
        if percent_aligned < float(min_percent_aligned):
            continue
   
        if gi_only:
            yield parse_gi_from_subject_id(blast_result[subject_index])
        elif add_sampleid:
            blast_result.append(add_sampleid)
            yield blast_result
        else:
            yield blast_result
        
def add_taxonomy_to_blast_results(lines,ncbi_taxonomy_tree,gi_to_taxid_mapping):
    """Output an annotated BLAST result with extra fields for taxid, and taxonomic lineage
    lines -- blast of a similarity search result 
       NOTE: the subject id is assumed to contain the gi in a format
        like this: gi|148686800|gb|EDL18747.1|
    
    header_map -- a dictionary mapping the names of columns to their index
    for example if column 0 is 'Query', header_map['Query']=0.
    
    ncbi_taxonomy_tree -- a PyCogent tree object representing the NCBI taxonomy tree
      made with PyCogents NCBITaxonomyTreeFromFiles, and NCBI's nodes.dmp and names.dmp
      files.
    
    gi_to_taxid_mapping -- a mapping between the relevant gis and taxids 
     
    """
    #The Subject field will always be in the same column
    #Look up ahead of time which column that is
    header_line,header_map,blast_result =\
      parse_similarity_search_lines(lines) 

    #Make a generator object that will yield parsed BLAST results 
    results_iterator =\
      yield_blast_results_with_taxonomy(blast_result,header_map,ncbi_taxonomy_tree,gi_to_taxid_mapping)
    
    #Update the header map and header line with new entries for 'gi','taxid',and 'taxonomy'
    last_field_idx = max(header_map.values())
    header_map["gi"] = last_field_idx +1
    header_map["taxid"] = last_field_idx +2
    header_map["taxonomy"] = last_field_idx +3
    new_header_line = header_line.strip()
    print "old header line:",header_line
    new_header_line = "\t".join([new_header_line,"gi","taxid","taxonomy"])
    print "new header line:",new_header_line
    return new_header_line,header_map,results_iterator


def yield_blast_results_with_taxonomy(results,header_map,ncbi_taxonomy_tree,gi_to_taxid_mapping):
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
    return lineage



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


def map_from_delimited_file(lines,key_field_idx=0,value_field_idx=1,delimiter="\t",\
    limit_keys_to=None):
    """ generates a mapping dictionary from two fields in a delimited file e.g. id1\tid2 
    lines -- an iterator of space delimited lines
    key_field_idx -- the index of the field that will be the keys to the dictionary
    value_field_idx -- the index of the field that will be the values in the dictionary
    limit_keys_to -- a set of key values.  keys not in the set
    """

    result  = {}
    for line in lines:
        fields = line.split(delimiter)
        curr_key = fields[key_field_idx]
        
        if limit_keys_to is not None:
            if curr_key not in limit_keys_to:
                continue

        curr_value = fields[value_field_idx]
        result[curr_key] = curr_value
    return result

def assign_sampleid_from_filename(input_fp):
    """Get a plausible default sampleid from the filename (removing extension)"""
    input_filepath,input_filename = path.split(input_fp)
    base_filename,extension = splitext(input_filename)
    return base_filename


def filter_and_merge_RAP_search_files(input_dir):
    #TODO: abstract step 1 here
    pass

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    
    #Load similarity results from file
    input_fp = opts.input_blast_results_file
    if not isdir(input_fp):
        if opts.verbose:
            print "Input is a single file: %s" %(input_fp)
        input_files = [input_fp]
    else:
        input_files = [f for f in listdir(input_fp)]
        if opts.verbose:
            print "%s is a directory. Scanning for files..." %(input_fp)
            print "Processing input files:","\n".join(listdir(input_fp))
    
    #Open output file now -- error out before processing if we can't open it
    output_dir = opts.output
    create_dir_return_code=create_dir(output_dir)
    if opts.verbose:
        print "Creating output directory:",output_dir
   
    #Load mapping between filenames and sampleids
    if opts.filename_to_sampleid_map_fp is not None:
        
        if opts.verbose:
            print "Loading filename to sampleid mapping"
        filename_to_sampleid_mapping_file=open(opts.filename_to_sampleid_map_fp,"U")
        filename_to_sampleid_mapping = map_from_delimited_file(filename_to_sampleid_mapping_file)
    else:
        filename_to_sampleid_mapping = None

    #Do Step 1.  Load a directory of RAP search results, with one file per sample
    filtered_RAP_file = open(output_dir+"/filtered_RAP_results.txt","w+")

    sampleids = set()
    for input_file in input_files:
    
        sampleid = None
        #Use a map file to get sampleids from filenames if one is supplied
        if filename_to_sampleid_mapping:
            input_filename,input_filepath = path.split(input_fp)
            sampleid = filename_to_sampleid_mapping.get(input_filename,None)

        #If no mapped sampleid was provided, generate a sampleid from filename
        if sampleid is None:
            sampleid = assign_sampleid_from_filename(input_file)

        if opts.verbose:
            print "File %s assigned to sampleid %s" %(input_file,sampleid)
            print "Adding sampleid to file"
    
        sampleids.add(sampleid)
    
        if opts.verbose:
            print "Filtering results by e-value and percent aligned" 
    
        #For the current file, filter results by e-value and percent id
        f = open("%s/%s" %(opts.input_blast_results_file,input_file),"U")
        filtered_RAP_result_iterator = filter_RAP_results(f,\
          opts.e_value,opts.percent_aligned,add_sampleid=sampleid,gi_only=False)
        for result in filtered_RAP_result_iterator:
            filtered_RAP_file.write("\t".join(result)+"\n")
        f.close()

    sampleids = list(sampleids)

    #Step 2. Filter to best hits only
    #For each filtered dataset (now of generally manageable size)
    #retain only the best RAP hit
    filtered_RAP_file.close()

    if opts.verbose:
        print "Done."
        print "The following sampleids were found: %s" %(sampleids)
        print "Filtering results by best BLAST hit"
    
    best_hits_file = open(output_dir+"/best_hits.txt","w")
    
    #We iterate over sampleids because we want to preserve
    #duplicate results across samples
    for sampleid in sampleids:
        if opts.verbose:
            print "Getting best hits for sampleid: %s" %sampleid
        
        filtered_results = open(output_dir+"/filtered_RAP_results.txt","U")
        
        header_line,header_map,best_hits =\
          best_RAP_results(filtered_results,\
          limit_to_sampleid = sampleid)

        filtered_results.close()
         
        #Write results to best hits file
        best_hits_file.write(header_line+"\n") 
        for best_hit in best_hits:
            best_hits_file.write("\t".join(best_hit)+"\n")
    
    
    best_hits_file.close() 
    
    #Reopen best hits file to read in gi numbers
    if opts.verbose: 
        print "Extracting GI numbers from best hits"
    best_hits_file = open(output_dir+"/best_hits.txt","U")
    relevant_gis = gis_from_blast_lines(best_hits_file)
    best_hits_file.close()

    if not relevant_gis:
        raise ValueError("No gi numbers could be extracted following filtering.  Try a less restrictive e-value or percent alignment?")
    #Step 3 Add taxonomy information
    if opts.verbose:
        print "Adding Taxonomy to BLAST results"
        print "Building NCBI Taxonomy Tree"
    
    ncbi_taxonomy_tree =\
      NcbiTaxonomyFromFiles(open(opts.ncbi_taxonomy_nodes_file,"U"),open(opts.ncbi_taxonomy_names_file,"U"))

    if opts.verbose:
        print "Building a GI to taxid map"

    #Because the gi->taxid files are very very large,
    #first create a set of relevant gis, then build a 
    #map for only those
    
    with open(opts.gi_to_taxid_file,"U") as f:
        gi_to_taxid_map = map_from_delimited_file(f,\
          limit_keys_to = relevant_gis)
    
    if opts.verbose:
        print "Done"
        print "Adding Taxonomy to BLAST results"

    best_hit_blast_lines = open(output_dir+"/best_hits.txt","U")
    
    header_line,header_map,annotated_blast_lines =\
      add_taxonomy_to_blast_results(best_hit_blast_lines,\
      ncbi_taxonomy_tree,gi_to_taxid_map)
   
    results_with_taxonomy_file = open(output_dir+"/best_hits_with_taxonomy","w")
    
    results_with_taxonomy_file.write(header_line+"\n")
    for result in annotated_blast_lines:
        results_with_taxonomy_file.write("\t".join(result)+"\n")
    results_with_taxonomy_file.close()
    best_hit_blast_lines.close()

    #Step 4 build an actual BIOM table
    
    with open(output_dir+"/best_hits_with_taxonomy","U") as f:
        result_biom_table = biom_table_from_blast_results(f)
    
    #Write BIOM output file to disk in JSON format
    with open(output_dir+"/otu_table.biom","w") as f:
        f.write(result_biom_table.getBiomFormatJsonString(generated_by=\
        "Viroscape version %s"%__version__))

    print "OTU table methods:",dir(result_biom_table)
    #Read in a file of blast results annotated with columns for taxid,sampleid
    #Numpy data matrix format: (observation,sample)

if __name__ == "__main__":
    main()
