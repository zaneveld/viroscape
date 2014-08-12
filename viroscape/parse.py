#!/usr/bin/env python
# File created on 2 Jul 2014
from __future__ import division


__author__ = "Jesse Zaneveld"
__copyright__ = "Vega Thurber Lab Copyright 2014"
__credits__ = ["Jesse Zaneveld"] 
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = ""
__email__ = ""
__status__ = "Development"


#Define a hardcoded header line for certain common similarity searches (e.g. RAP search -m 8).
RAP_HEADER_LINE = "# Fields:Query\tSubject\tidentity\taln-len\tmismatch\tgap-openings\tq.start\tq.end\ts.start\ts.end\tlog(e-value)\tbit-score\n" 

def parse_similarity_search_lines(lines,strict=False,default_header=RAP_HEADER_LINE):
    """Parse similarity search (e.g. BLAST or RAP search results) in tab-delimited (-m 8) format
    lines -- lines of results
    strict -- if True, error if no header line is found
    default_header -- if a string is provided, use it as a default headerline if 
      no header is present in file. Header lines must start with the string '# Fields:'
      and occur before any non-comment lines.
    
    Comment lines starting with '#' will be ignored
    Returns the header line (from the data or supplied default), a dict mapping header fields
      to their index in the field list (e.g. e-value is column 7, {'e-value':7}, and a generator
      that yields lists of fields for all of the data lines.
    """
    #First parse comment lines to find a header 
    header_map = None 
    header_line = None 
    for i,line in enumerate(lines): 

        #Header line must before any data lines 
        if line.startswith("# Fields:"): 
            header_line = line 
            header_map = make_header_map(header_line) 
            break 
        elif line.startswith("#"): 
            #Skip other comment lines 
            continue 
        elif header_map is None and not strict and default_header is not None: 
            #If we have reached a non-comment line without finding a header
            #there must not be one in this file.  
             
            #Use the default header line provided by the user
            header_line = default_header 
            header_map = make_header_map(header_line) 
            break  
        elif header_map is None: 
            #We must be in a data line, but never 
            #founda  header. 
            #Either the user set strict = True or no default header was supplied.
            raise ValueError("Never found a header line starting with '# Fields:'") 
    
    #Next parse actual data lines to give a generator of data fields
    data_fields = parse_delimited_data_lines(lines) 
    
    return header_line,header_map,data_fields 
 
def parse_delimited_data_lines(data_lines,delimiter="\t"):
    """Parse blast data (non-header) lines, yielding lists of fields
    
    Blank lines and comment lines starting with a '#' sign will be skipped
    """
    for line in data_lines: 
         
        if line.startswith("#"): 
            continue 
        if not line.strip(): 
            continue 
         
        fields = line.strip().split(delimiter) 
        yield fields 
 
def parse_gi_from_subject_id(subject_id): 
    """Parse a numerical gi from a subject id of the form: gi|148686800|gb|EDL18747.1|""" 
    fields = subject_id.split("|") 
    for i in xrange(0,len(fields),2): 
        paired_data = fields[i:i+2] 
        database_name,identifier = paired_data 
        if database_name == 'gi': 
            return identifier 
    
def make_header_map(header_line):
    """Yield a dict mapping fields to headers"""
    # Fields: Query Subject identity    aln-len mismatch    gap-openings    q.start q.end   s.start s.end   log(e-value)    bit-score 
    if not header_line.startswith("# Fields:"):
        raise ValueError("Results header line must look like this:\n # Fields: Query Subject identity    aln-len mismatch    gap-openings    q.start q.end   s.start s.end   log(e-value)    bit-score")
    header_line = header_line.lstrip("# Fields:").strip()
    header_map = {key:idx for idx,key in enumerate(header_line.split("\t"))}
    return header_map


