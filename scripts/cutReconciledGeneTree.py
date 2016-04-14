#!/usr/bin/env python
__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__credits__ = ["Mehmood Alam Khan"]
'''
Created on Nov 10, 2014

@author: malagori
'''

from optparse import OptionParser
from dendropy import Tree
import sys
import os
import subprocess
import random
import tempfile
import csv
import glob


def checkExe(exePath):
    return os.path.isfile(exePath) #and os.access(exePath, os.X_OK)
    
def Where(program):
    '''
    input: name of executable
    output: path to executable
    '''
    fpath, fname = os.path.split(program)
    if fpath:
        if checkExe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            pathToProgram = os.path.join(path, program)
            if checkExe(pathToProgram):
                return pathToProgram
    return None

def allLeafLabels( myTree):
    '''
    returns the list containing all the leaf nodes' label of this tree and the internal nodes without the root label
    '''
    leafLabels=[]
    internalLabelsExRoot=[]
    for i in myTree.leaf_nodes():
        leafLabels.append(i.get_node_str().replace("'", ""))
    for i in myTree.internal_nodes():
        if i.level() != 0:
            internalLabelsExRoot.append(i.label)
    return leafLabels, internalLabelsExRoot

def onlyLeafLabels( myTree):
    '''
        returns the list containing labels of leaf nodes for this tree
    '''
    leafLabels=[]
    for i in myTree.leaf_nodes():
        leafLabels.append(i.get_node_str().replace("'", ""))
    return leafLabels

def readTreeFromFile( treePath):
    '''
    input: path to the file containing newick tree
    return Tree object 
    '''
    print treePath
    #myTree= Tree()
    myTree= Tree.get_from_path(treePath, 'newick', annotations_as_nhx=True, extract_comment_metadata=True , suppress_annotations=False)
    return myTree

def readTreeFromString( treeString):
    '''
    input: string containing newick tree
    return Tree object 
    '''
    #myTree= Tree()
    myTree= Tree.get_from_string( treeString, 'newick', annotations_as_nhx=True, extract_comment_metadata=True , suppress_annotations=False)
    
    return myTree
 
def readLeafmapMappingFile(infile):
    '''
    infile format: primeGeneID <tab> newGeneID <tab> newSpeciesID
    '''
    primeGeneIdMappingDic= {}
    with open(infile, 'r') as rf:
        for line in rf:
            line= line.rstrip()
            tok=[]
            tok= line.split()
            primeGeneIdMappingDic[tok[0]]=str(tok[2]+'_'+tok[1])
    return primeGeneIdMappingDic
        
def readGuest2HostMappingFile(infile):
    '''
    infile format: Guest vertex ID: <tab> Host vertex/arc ID:
    '''
    guest2HostMappingDic={}
    with open(infile, 'r') as rf:
        next(rf) # skip the first line in file
        for line in rf:
            line= line.rstrip()
            tok=[]
            tok=line.split()
            guest2HostMappingDic[tok[0]]= tok[1]
    return guest2HostMappingDic
            
def readRecGtreeAndCut( gTree, workDir):
    '''
        This function take a single reconciliation and cut at speciation nodes
    '''
    orthoFamilies       = {} # key= root label of subtree; value: subtree
    speciesOrthoCuts    = {}
    recGTree= readTreeFromFile(gTree)
    mappingFileNumber=gTree.split('/')[-1].split('.')[0]
    mappingFileDir=os.path.dirname(os.path.abspath(gTree))
    mFile=mappingFileDir+'/'+mappingFileNumber+'.pruned.leafmap.mapping'
    guest2HostMappingFileName= mappingFileDir+'/'+str(mappingFileNumber)+'.guest2hostmapping'
    primeGeneIdMappingDic=readLeafmapMappingFile(mFile)
    cutAtSpeciationVertices( orthoFamilies, guest2HostMappingFileName, recGTree, speciesOrthoCuts, primeGeneIdMappingDic)

    try:
        tok=[]
        tok= gTree.split('/')
        orthoHandle = csv.writer(open(workDir+'/'+tok[-1]+".sallfamilies",  'wb'))
        for key, value in orthoFamilies.items():
            orthoHandle.writerow([key, value, speciesOrthoCuts[key]])
    except IOError, e:
        print ("Class: Wrapper, Function: readRecGtreeAndCut(): %s " % e)



def cutAtSpeciationVertices( orthoFamilies, guest2HostMappingFileName, recGTree, speciesOrthoCuts, primeGeneIdMappingDic):
    '''
        This function cuts the reconciled gene tree at the speciation vertices.
        input:     Reconciled Gene Tree, path to the file containing species tree
        output:    families. soft orthologes and soft paraloges
    '''
    
    guest2HostMappingDic=readGuest2HostMappingFile(guest2HostMappingFileName)
    
    for i in recGTree.internal_nodes():
        leafLabels      =[]
        idVtypeDiscPT   =[]
        vType           =[]
        idVtypeDiscPT= i._annotations.values_as_dict()['PRIME ID'].split()
        primeID=idVtypeDiscPT[0]
        cutNode=guest2HostMappingDic[primeID]
        vType= idVtypeDiscPT[1].split("=")
        if  vType[1] == 'Speciation':
            # orthologous families
            t= Tree(i)
            leafLabels= onlyLeafLabels(t)
            newLeafLabels=[primeGeneIdMappingDic[k] for k in leafLabels]
                
            family=""
            family= str(sorted(newLeafLabels))
            
            
            if cutNode not in speciesOrthoCuts.keys():
                speciesOrthoCuts[family]=cutNode
            
            if family not in orthoFamilies.keys():
                orthoFamilies[family]=1
            else:
                orthoFamilies[family] +=1
                
     
                    
def Cutter( gTree, workDir):
    '''
        This function control the work flow of this simulator.
    '''
    readRecGtreeAndCut( gTree, workDir)
    
    
def main(argv):
    #global leafLabelStree, internalExRootSpeceLabels
    
    '''
    This program takes gene tree generated by phyloGenData program and cut the tree on speciations events
    '''
   
        
    usage = "usage: ./cutReconciledGeneTree [options]"
    
    description="""DESCRIPTION:     This program takes gene tree generated by phyloGenData program and cut the tree on speciations events"""

    if len(argv) == 0:
        print usage
        sys.exit()
        
    parser = OptionParser(usage=usage, description=description )
    
    parser.add_option('-d','--genesDir', 
                      action="store", type="string", dest="genesDir", 
                      help='Specify path to the directory containing gene trees.'+ 
                      ' Note: this option is used when the input gene trees are more than one.')
    parser.add_option('-p','--geneTreeExtension', 
                      action="store", type="string", dest="geneTreeExtension", 
                      help='Specify the extension of gene tree name. e.g: .newick or .tree etc'+
                      ' Note: This option is used with -d option.')
    parser.add_option('-G','--geneTree', 
                      action="store", type="string", dest="geneTree", 
                      help='Specify path to the gene tree file.')
    parser.add_option('-D', '--outDirectory',
                      action="store", type="string", dest="outDir",
                      help='Specify path to directory to store intermediate results. If not specified intermediate results will be stored in a temporary directory which is deleted after program execution (default=None)', default=None)

    (options, args) = parser.parse_args()
    
    
    
    if options.outDir == None:
        options.outDir=tempfile.mkdtemp()
    elif os.path.exists(options.outDir) == False:
        os.mkdir(options.outDir)
    print "Output directory: %s"  % (options.outDir)

    if options.genesDir != None:
        print "Gene tree Directory: %s"       % (options.geneTree)
        if options.geneTreeExtension == None:
            print "Input Error: Please provide Gene Tree Extension i.e. use -p option"
            sys.exit()
        else:
            allReconciledGeneTrees= glob.glob(options.genesDir+"/*."+options.geneTreeExtension)
            for i in allReconciledGeneTrees:
                print "Gene Tree: %s" %(str(i))
                Cutter( i, options.outDir)
    else:
        if checkExe(options.geneTree) == True:
            print "Gene tree: %s"       % (options.geneTree)
        else:
            print "Error: Gene tree file not found"
            sys.exit()
        Cutter( options.geneTree, options.outDir)
    

if __name__ == '__main__':
    main(sys.argv[1:])
