#!/usr/bin/env python
'''
Created on April 07, 2014

@author: Mehmood Alam Khan Malagori
@email:   malagori@kth.se
'''
__version__= "0.1dev"
__license__ = "GPLv3"

import argparse
import sys
from dendropy import Tree

def changeSpeciesTreeLabels(stree):
    '''
        this function change the labels of species tree to new names
    '''
    myStree= Tree.get_from_path(stree, 'newick', annotations_as_nhx=True, extract_comment_metadata=True , suppress_annotations=False)
    myStree.print_plot()
    k=0
    with open(stree+'.labels', 'w') as wf:
        for n in myStree.leaf_nodes():
            #wf.write(n.taxon.label +'\t'+ 'S'+str(k+1) +'\n')
            wf.write(n.taxon.label +'\t'+ str(k+1) +'\n')
            #n.taxon.label= 'S'+str(k+1)
            n.taxon.label= str(k+1)
            k=k+1
    myStree.print_plot()
    with open(stree+'.newNewick', 'w') as wf:
        st=myStree.as_string('newick')
        wf.write(st)
        
def main(argv):
    
    # Take care of input
    parser = argparse.ArgumentParser(description="Parse input arguments and print output.")
    
    parser.add_argument('-s', metavar='sTree',type=str, help='Specify path to the species tree file to be converted into new species labels ', default= None)
    
    args = parser.parse_args()
    stree= args.s
    
    changeSpeciesTreeLabels(stree)


if __name__== "__main__":
    main(sys.argv[1:])