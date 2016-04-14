#!/usr/bin/env python
__author__ = "Mehmood Alam Khan"
__email__  = "malagori@kth.se"
__version__= "0.9"
__credits__ = ["Mehmood Alam Khan"]
'''
Created on Sep 24, 2014

@author: malagori
'''
from optparse import OptionParser
#from dendropy import Tree
import csv
import sys
import os
import random
from ete2 import Tree, TreeStyle, NodeStyle, faces, AttrFace, BarChartFace

def layout(node):

    if node.is_leaf():
        # Add node name to laef nodes
        N = AttrFace("name", fsize=20, fgcolor="black")
        faces.add_face_to_node(N, node, 0)
    elif "pie_data" in node.features:
        # Creates a pie chart face whose size is proportional to node's
        percents=node.pie_data
        rcolor=[]
        for i in xrange(0, len(percents)):
            r= lambda: random.randint(0,255)
            rcolor.append( ('#%02X%02X%02X' % (r(),r(),r())))
        C = BarChartFace(percents, width=WIDTH, height=HEIGHT, colors=rcolor, min_value=0, max_value=100)
        # Let's make the Bar Chart transparent 
        #C.opacity = 0.3
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, column=0, position="branch-top")
        nameFace= AttrFace("name", fsize=20)
        faces.add_face_to_node(nameFace, node, column=0, position="branch-bottom")
def readTreeFromFile(treePath):
    '''
    input: path to the file containing newick tree
    return Tree object 
    '''

    myTree= Tree(treePath, format=1)
    myTree.convert_to_ultrametric(myTree.__len__(), strategy='balanced')
    return myTree

def getFamiliesStatisticsForEachNode(pathToSfamilies, bootValue):
    """
    This function takes the path to the sfamilies file and return the statistics over each node of species tree.
    output: snodes dic: key=nodes_name, values= list of families support values for each speceis node's cut
    """
    snodes={}
    sumAllSupport=0
    with open(pathToSfamilies, 'r') as rf:
        sFamReader= csv.reader(rf)
        
        for i in sFamReader:
            if i[2] not in snodes.keys():
                snodes[i[2]]=list()
                snodes[i[2]].append(int(i[1])) # storing the support for each family not string family
               
                family=i[0].replace("[",'')
                family=family.replace("]",'')
                family=family.replace(",", ' ')
                family=family.replace("'", '')
                family=family.split()
                familySize= len(family)
            else:
                snodes[i[2]].append(int(i[1])) # storing the support for each family not string family
               
                family=i[0].replace("[",'')
                family=family.replace("]",'')
                family=family.replace(",", ' ')
                family=family.replace("'", '')
                family=family.split()
                familySize= len(family)
    # converting support values to range 0 to 100
    
    oldMax= int(bootValue)
    oldMin= 0
    newMax= 100
    newMin= 0
    
    oldRange = (oldMax - oldMin)  
    newRange = (newMax - newMin)  
    
    for key, values in snodes.iteritems():
        newValueList=[]
        for i in values:
            newValue = (((i - oldMin) * newRange) / oldRange) + newMin
            newValueList.append(newValue)
            sumAllSupport += newValue
        snodes[key]=newValueList
       
    #for key, values in snodes.iteritems():    
    #    newVal=[]
    #    for i in values:
    #        v= 1/float(sumAllSupport) * 100
    #        newVal.append(v)
    #    snodes[key]= newVal

    return snodes
            
    
def visualizeTree(sTreePath, pathToSfamilies, bootValue, width, height):
    # Random tree
    stree = Tree()
    stree = readTreeFromFile(sTreePath)
   
    snodesStatDic={}
    snodesStatDic= getFamiliesStatisticsForEachNode(pathToSfamilies, bootValue)
    #print snodesStatDic
    # Some random features in all nodes
    for n in stree.traverse():
        if n.name in snodesStatDic.keys():
            total= reduce(lambda x,y: x+y, snodesStatDic[n.name])
            #norm= [(x*100)/total for x in snodesStatDic[n.name]]
            norm= [x for x in snodesStatDic[n.name]]
            n.add_features(pie_data=norm)
    # Create an empty TreeStyle
    ts = TreeStyle()

    # Set our custom layout function
    ts.layout_fn=layout

    # Draw a tree 
    ts.mode = "r"
    
    #ts.force_topology= False
    ts.complete_branch_lines_when_necessary= True
    # We will add node names manually
    ts.show_leaf_name = False
    # Show branch data
    #ts.show_branch_length = True
    #ts.show_branch_support = True
    

    return stree, ts


def checkExe(exePath):
    return os.path.isfile(exePath)

def main(argv):
    global WIDTH, HEIGHT
    
    '''
    This program takes the species tree and sfamilies output file from GFI program.
    It outputs the results in a visualized form.
    '''
   
        
    usage = "usage: ./visualizeResults [options]"
    
    description="""DESCRIPTION:     This program takes the species tree and sfamilies output file from GFI program.
    It outputs the results in a visualized form."""

    if len(argv) == 0:
        print usage
        sys.exit()
        
    parser = OptionParser(usage=usage, description=description )
    
    parser.add_option('-S','--speciesTree', 
                      action="store", type="string", dest="speciesStree", 
                      help='Specify path to the species tree file.')
    parser.add_option('-f','--sfamilies', 
                      action="store", type="string", dest="sfamilies", 
                      help='Specify path to the sfamilies file.')
    parser.add_option('-b','--bootValue', 
                      action="store", type="string", dest="bootValue", 
                      help='Specify the bootstrap value used for this sfamilies file.')
    parser.add_option('-w','--width', 
                      action="store", type="int", dest="width", 
                      help='Specify the width of barchart, default:200.', default=200)
    parser.add_option('-H','--height', 
                      action="store", type="int", dest="height", 
                      help='Specify the height of barchart, default:200.', default=200)
    parser.add_option('-o','--outputFile', 
                      action="store", type="string", dest="outputFile", 
                      help='Specify path to the outputFile file.')
    
    
    (options, args) = parser.parse_args()
    

    if checkExe(options.speciesStree) == True:
        print "Species tree: %s"    % (options.speciesStree)
    else:
        print "Error: Species tree file not found"
        sys.exit()
    if options.speciesStree == None or options.bootValue == None or options.sfamilies == None or options.outputFile == None:
        print "Too few arguments!"
        sys.exit()
    WIDTH= options.width
    HEIGHT= options.height
    stree, myTreeStyle= visualizeTree(options.speciesStree, options.sfamilies, options.bootValue, options.width, options.height)
    #print dir(stree.show)
    #stree.set_style(myTreeStyle)
    #stree.show(layout)
    #stree.show(tree_sytle=myTreeStyle)
    stree.render(options.outputFile, w=500, h=500, dpi=300, tree_style=myTreeStyle)

if __name__ == '__main__':
    main(sys.argv[1:])
