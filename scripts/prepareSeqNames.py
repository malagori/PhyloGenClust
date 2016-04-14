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

def fixSequenceNames(seqFile,directory, gsMap, gTree):
    '''
     This function create sequence labels such that its consistence with Notung preorder functionality
    '''
    tok=[]
    seqLength=0
    #tok=seqFile.split('.')
    #newFile= str(tok[0]+tok[1])
    newFile= seqFile+".prefixed"
    speciesNames={}
    genesMapping={}
    speciesMapping={}
    S=0
    G=0
    #uniqueSpeciesList=[]
    with open(directory+"/"+newFile, 'w') as wf, open(directory+"/"+gsMap+".mapping", 'w') as mwf:
        with open(str(directory+"/"+gsMap), 'r') as gsf:
            for line in gsf:
                tok=[]
                tok=line.rstrip().split('\t')
                #speciesTok=[]
#                speciesTok=tok[1].rstrip().split('_') # split species name delimited by _
#                speciesNameList=[x.capitalize() for x in speciesTok]
#                speciesName= "".join(speciesNameList)
#                speciesName= speciesName.strip()
                
#                if tok[-1] not in speciesMapping.keys():
#                    speciesMapping[tok[-1]] = S+1
#                    if tok[0] not in genesMapping.keys():
#                        genesMapping[tok[0]] = G+1 
#                        speciesNames[G+1]= S+1
#                        G= G+1
#                        S= S+1
#                    else:
#                        speciesNames[genesMapping[tok[0]]]= S+1
#                        S= S+1
#                else:
#                    if tok[0] not in genesMapping.keys():
#                        genesMapping[tok[0]] = G+1 
#                        speciesNames[G+1]= speciesMapping[tok[-1]]
#                        G= G+1
#                    else:
#                        speciesNames[genesMapping[tok[0]]]= speciesMapping[tok[-1]]

                if tok[0] not in genesMapping.keys():
                    genesMapping[tok[0]] = G+1
                    speciesNames[G+1]= tok[-1]  
                    G = G+1
                else:
                    speciesNames[genesMapping[tok[0]]]= tok[-1]
                # gene name, new gene name, species name
                mwf.write(tok[0]+ "\t"+ str(genesMapping[tok[0]])+ "\t"+tok[-1]+"\n" )                 
                #speciesNames[tok[0]]=tok[-1]
            
        flag=True
#        uniqueSpecies= set(speciesNames.values())
#        uniqueSpeciesList= list(uniqueSpecies)
        with open(str(directory+'/'+seqFile), 'r') as seqf: 
            for i in seqf:
                if flag == False:
                    tok=[]
                    tok=i.strip().split(" ")
                    geneName=tok[0].strip()
                    newGeneName= str(speciesNames[genesMapping[geneName]])+'_'+str(genesMapping[geneName])
                    newRecord= str('{:10s} {:'+str(seqLength)+'}').format(newGeneName, tok[-1])

                    wf.write(newRecord+'\n')
                else:
                    flag = False
                    tok=[]
                    tok=i.strip().split(" ")
                    seqLength=tok[-1]
                    wf.write(i)
    addSpeciesPrefixToGTree(directory, gTree, speciesNames, genesMapping)

def readTreeFromFile( treePath):
    '''
    input: path to the file containing newick tree
    return Tree object 
    '''
    myTree= Tree.get_from_path(treePath, 'newick', annotations_as_nhx=True, extract_comment_metadata=True , suppress_annotations=False)
    return myTree

def addSpeciesPrefixToGTree(directory, gTree, speciesNames, genesMapping):
    """
    This function will add species name at the beginning of gene name to all the gene trees in the specified directory.
    e.g 1_10
        1_2
        ...
    """
    
    gT= readTreeFromFile(str(directory+"/"+gTree))
    for n in gT.leaf_nodes():
        tok = []
        tok = n.taxon.label.split()
        #n.taxon.label= speciesNames[tok[0]+'_'+tok[1]]+'_'+tok[0]+'_'+tok[1]
        n.taxon.label= str(speciesNames[genesMapping[tok[0]+'_'+tok[1]]])+'_'+ str(genesMapping[tok[0]+'_'+tok[1]])
    #gT.print_plot()
    
    #gT.write_to_path(directory+"/"+gTree+".prefixed", 'newick', suppress_internal_node_labels=False,annotations_as_nhx=True, extract_comment_metadata=True , suppress_annotations=False)
    
    with open(str(directory+"/"+gTree+".prefixed"), 'w') as wf:
        st=gT.as_string('newick').replace("'", "")  # this is really important to remove ' form node names otherwise notung will not read
        wf.write(st)

    
def main(argv):
    
    # Take care of input
    parser = argparse.ArgumentParser(description="Parse input arguments and print output. works only for phylip format data.")
    
    parser.add_argument('-d', metavar='dir',type=str, help='Specify path to the directory containing gene trees, sequences, gs map files in same directory ', default= None)
    parser.add_argument('-n', metavar='nGeneTrees',type=int, help='Specify the number of gene trees in the given directory ', default= None)
    #parser.add_argument('-gd', metavar='genesdirectory',type=str, help='Specify path to the directory containing gene trees in *.relaxed.tree format, otherwise you need to change in the code ', default= None)
    #parser.add_argument('-l', metavar='seqLength',type=int, help='Specify the length of sequences ', default= None)
    #parser.add_argument('-s', metavar='sTree',type=str, help='Specify path to the species tree file ', default= None)
    
    args = parser.parse_args()
    
    #seqFile = args.d
    #gsMap   = args.gs
    #seqLength= args.l
    #directory= args.gd
    #stree= args.s
    
    directory= args.d
    nGTrees= args.n
    
    for i in xrange(1, nGTrees+1):
        #if i!=302:
        fixSequenceNames(str(i)+".phylip",directory, str(i)+".pruned.leafmap",str(i)+".relaxed.tree")

if __name__== "__main__":
    main(sys.argv[1:])