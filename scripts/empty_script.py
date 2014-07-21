#!/usr/bin/env python
# File created on 2 Jul 2014
from __future__ import division


#of credits
__author__ = ""
__copyright__ = "Vega Thurber Lab Copyright 2014"
__credits__ = ["",""] 
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = ""
__email__ = ""
__status__ = "Development"
 

from cogent.util.option_parsing import parse_command_line_parameters, make_option
from os import path

script_info = {}
script_info['brief_description'] = "This script ... <TODO: add a one sentence explanation of what the script does here>"
script_info['script_description'] = "This is an empty template script useful as a starting point for adding a nice interface to your python scripts. PyCogent is the only dependency.  <TODO:remove thsi filler text and  add a longer explanation of the script here>"
script_info['script_usage'] = [("","< examples of script usage go here>. Predict metagenome from count file and OTU table.","%prog -i otus.biom -c KEGG_acepic__predict_traits_97.biom.gz -o predicted_metagenomes.biom"),
                               ("","<You can have more than one example>. Change output format to plain tab-delimited:","%prog -f -i otus.biom -c KEGG_acepic_predict_traits_97.biom.gz -o predicted_metagenomes.tab")]
script_info['output_description']= "<Describe output data and format here> Output is a table of function counts (e.g. KEGG KOs) by sample ids."
script_info['required_options'] = [
 make_option('-i','--input_file',type="existing_filepath",help='<TODO: describe input data here>. Example: the input trait counts on per otu basis in biom format (can be gzipped)'),
 make_option('-o','--output',type="new_filepath",help='<TODO: describe output here> The output filepath for results')
]
script_info['optional_options']=[
 make_option('-x','--example_optional_option',action="store",type="int",default=0,help='the output filepath for the scaled metagenome')
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

if __name__ == "__main__":
    main()
