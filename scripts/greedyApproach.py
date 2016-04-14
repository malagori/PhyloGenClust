#!/usr/bin/env python
'''
Created on Feb 3, 2015

@author: malagori
'''

import sys
import os
import glob
import csv
import operator
from dendropy import Tree
from optparse import OptionParser


def returnRootOfTree( infile, filePrefix, ext):
    '''
    input: path to the file containing newick tree
    return root of the Tree 
    '''
    
    directory=os.path.dirname(os.path.realpath(infile))
    treePath= directory+'/'+filePrefix+'.'+ ext
    rootNode=''
    myTree= Tree.get_from_path(treePath, 'newick', annotations_as_nhx=True, extract_comment_metadata=True , suppress_annotations=False)
    for i in myTree.internal_nodes():
        if i.level() == 0:
            rootNode=i.get_node_str()
            break
    return rootNode

def removeExtraCharaters(dirtyString):
    j=dirtyString.replace("[",'')
    j=j.replace("]",'')
    j=j.replace(",", ' ')
    j=j.replace("'", '')
    j=j.split()
    
    return j
    

def mergeGreedly(lowerThreshold, upperThreshold, sallfamDir, outDir, IGNORNEGENES, ext):
    '''
    Apply threshold to families
    returns sub families for each big family passes the threshold criteria
    '''
    # Pick families which are above the given threshold
    file_count=0
    allOrthoFilesList= glob.glob(sallfamDir+"/*.sallfamilies")

    # write all orthologous families above the given threshold
    #with open(outDir+"/results."+str(lowerThreshold)+"_"+str(upperThreshold)+".greedy", 'w') as fwf, open (outDir+ '/root_list.csv', 'w') as rootf:
    with  open (outDir+ '/root_list.csv', 'w') as rootf:    
        #allFamWriter= csv.writer(fwf)
        #allFamWriter.writerow(["ClusterId", "PredictedFamily", "SizeOfPredictedFamily", "Counts", "CutId"])
        
        for orthoFile in allOrthoFilesList:
                        
            twoDdic={}
            fname= orthoFile.split("/")[-1]
            f= fname.split('.')
            
#            rot=''
#            rot= returnRootOfTree(orthoFile, f[0], ext)
#            rootf.write("familyID: %s, rootName: %s \n" % (sallfamDir+"."+f[0], rot))
#            
            #with open (orthoFile, 'rb') as rf, open (truthFile, 'rb') as tf, open(outDir+"/"+fname+"."+str(lowerThreshold)+"_"+str(upperThreshold)+".greedy", 'w') as wf:
            with open (orthoFile, 'rb') as rf, open(outDir+"/"+fname+"."+str(lowerThreshold)+"_"+str(upperThreshold)+".union", 'w') as wf, open(outDir+"/"+fname+"."+str(lowerThreshold)+"_"+str(upperThreshold)+".intersect", 'w') as iwf:
            
                unionFamiliesWriter = csv.writer(wf)
                interFamiliesWriter = csv.writer(iwf)
                orthoReader= csv.reader(rf)
                #truthReader= csv.reader(tf)
                
                unionFamiliesList    = []
                intersectFamiliesList= []
                greedyUnionAtCutNodes={}
                greedyIntersectAtCutNodes={}
                
                for i in orthoReader:
                    twoDdic.setdefault(i[2],{})[i[0]]=i[1]  # {'SpeciesCutNode': {'PredictedFamily': bootNum, 'PredictedFamily': bootNum, ...}}
                
                for key, value in twoDdic.items():
                    familyBootnum = sorted(value.items(), key=operator.itemgetter(1))
                    familyBootnum.reverse()
                    unionFamiliesList    = []
                    intersectFamiliesList= []
                    
                    
                    
                    for i in xrange(0, len(familyBootnum)):
                        if int(familyBootnum[i][1]) > int(lowerThreshold) and int(familyBootnum[i][1]) < int(upperThreshold):
                                
                            j=removeExtraCharaters(familyBootnum[i][0])
                            for m in xrange(i+1, len(familyBootnum)):
 
                                #if int(familyBootnum[m][1]) > int(lowerThreshold) and int(familyBootnum[m][1]) < int(upperThreshold):
                                if int(familyBootnum[m][1]) < int(upperThreshold):  
                                    
                                    #print "here"
                                    k= removeExtraCharaters(familyBootnum[m][0])
                                    unionFamiliesList.append(removeExtraCharaters(familyBootnum[i][0]))
                                    intersectFamiliesList.append(removeExtraCharaters(familyBootnum[i][0]))
                                    
                                    diffGenes= list(set(k).symmetric_difference(set(j)))
                                    
                                    if len(diffGenes) <= round(IGNORNEGENES*len(k)): # 30% of the members will be ignored
                                        unionFamiliesList.append(sorted(list(set(k).union(set(j)))))
                                        intersectFamiliesList.append(sorted(list(set(k).intersection(set(j)))))
                    greedyUnionAtCutNodes[key]= unionFamiliesList
                    greedyIntersectAtCutNodes[key]= intersectFamiliesList
                uniqueUion=[]
                uniqueInter=[]
                for l, v in greedyIntersectAtCutNodes.iteritems():
                    for x in v:
                        if x not in uniqueInter:
                            #print x, l
                            uniqueInter.append(x)
                            interFamiliesWriter.writerow([x,l])
                            
                for l, v in greedyUnionAtCutNodes.iteritems():
                    for x in v:
                        if x not in uniqueUion:
                            
                            uniqueUion.append(x)
                            #print x, l
                            unionFamiliesWriter.writerow([x,l])
#            file_count +=1
#            if (file_count == 1):
#                break
        
def main(argv):

    '''
    Main program
    '''
    usage = "usage: ./greedyApproach.py [options]"
    if len(argv) == 0:
        print usage
        sys.exit()
    parser = OptionParser(usage=usage)
    parser.add_option('-f','--sallfamiliesDir', 
                      action="store", type="string", dest="sallfamiliesDir", 
                      help='Specify path to the directory containing all the *.allfamilies file.')
    
    parser.add_option('-l','--lowerThreshold', 
                      action="store", type="float", dest="lowerThreshold", 
                      help='Provide the threshold value range[0, 1]. default: 0.3', default=0.3)
    parser.add_option('-u','--upperThreshold', 
                      action="store", type="float", dest="upperThreshold", 
                      help='Provide the threshold value range[0, 1]. default: 0.5', default=0.5)
    parser.add_option('-i','--ignoreNumGene', 
                      action="store", type="float", dest="ignoreNumGene", 
                      help='Specify the number of conflicting genes to be ignored while merging two families.default=20% of family size', default=0.2)
    
    parser.add_option('-e','--extension', 
                      action="store", type="string", dest="extension", 
                      help='Provide the extension of species Tree.')
    parser.add_option('-b','--bootValue', 
                      action="store", type="int", dest="bootValue", 
                      help='Specify the bootstrap value used for this sfamilies file.')
    parser.add_option('-D', '--outDirectory',
                      action="store", type="string", dest="outDir",
                      help='Specify path to directory to results.', default=None)

    
    (options, args) = parser.parse_args()
    
    if os.path.exists(options.outDir) == False:
        os.mkdir(options.outDir)
        
    
    if options.sallfamiliesDir == None or options.bootValue == None or options.lowerThreshold == None or options.upperThreshold == None or options.outDir == None: #options.nFamilies == None:
        print "Too less arguments"
        sys.exit()
    print "Writing files to output directory: %s"  % (options.outDir)
    lt= options.lowerThreshold * options.bootValue
    up= options.upperThreshold * options.bootValue
    mergeGreedly(lt, up,  options.sallfamiliesDir, options.outDir, options.ignoreNumGene , options.extension)    
    
    
if __name__ == '__main__':
    main(sys.argv[1:])