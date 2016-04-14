'''
Created on Sep 19, 2013

@author: Mehmood Alam Khan Malagori
@email:   malagori@kth.se
'''
__version__= "0.1dev"
__license__ = "GPLv3"


import subprocess
import os
import sys
import random as rn
import csv

from dendropy import Tree
from Bio.Emboss.Applications import FSeqBootCommandline
from Bio.Phylo.Applications import FastTreeCommandline
from cStringIO import StringIO


class Wrapper(object):
    '''
    This class contains wrappers for all the tools used in this project
    '''


    def __init__(self, workDir, resultsFile, inFile, coreId, seqType, seedNum, bootNum, method, interLeaved):
        '''
        Data Fields:
        work_dir   = temproray directory
        inFile     = data file,
        id         = core id (int), 
        seqType    = d (dna); p (protein); r (rna), 
        bootNum    = number of replicates,
        seedNum    = Random number seed between 1 and 32767
        method     =b (Bootstrap) Default, 
                    j (Jackknife)
                    c (Permute species for each character)
                    o (Permute character order)
                    s (Permute within species)
                    r (Rewrite data)),
        interLeaved=True if sequence data is interleaved otherwise False
        '''
        self.work_dir       = workDir
        self.resultsFile    = resultsFile
        self.inFile         = inFile 
        self.coreId         = coreId
        self.seqType        = seqType
        self.seedNum        = seedNum
        self.bootNum        = bootNum
        self.method         = method
        self.n              = bootNum
        self.nt             = False
        self.interLeaved    = interLeaved
        self.outFile        = "bootstrap_"+str(coreId)+".out"
        self.newSpeciesTree = Tree()
        self.leafLabelStree = []
        self.internalExRootSpeceLabels= []
        
        if self.seqType in ['r', 'd']:
            self.nt = True
    
    def getSpeciesTree(self):
        return self.newSpeciesTree
    
    def setSpeceistree(self, stree):
        self.newSpeciesTree= stree

    def writeSpeciesTree(self):
        self.newSpeciesTree.write_to_path(self.resultsFile+'.speciestree.'+str(self.coreId), 'newick', suppress_internal_node_labels=False,annotations_as_nhx=True, extract_comment_metadata=True , suppress_annotations=False)
    
    def allLeafLabels(self, myTree):
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
    def onlyLeafLabels(self, myTree):
        '''
            returns the list containing labels of leaf nodes for this tree
        '''
        leafLabels=[]
        for i in myTree.leaf_nodes():
            leafLabels.append(i.get_node_str().replace("'", ""))
        return leafLabels
        
    def readTreeFromFile(self, treePath):
        '''
        input: path to the file containing newick tree
        return Tree object 
        '''
        print treePath
        myTree= Tree()
        myTree= Tree.get_from_path(treePath, 'newick', suppress_edge_lengths=False, annotations_as_nhx=True, extract_comment_metadata=True , suppress_annotations=False)
        return myTree

    def readTreeFromString(self, treeString):
        '''
        input: string containing newick tree
        return Tree object 
        '''
        myTree= Tree()
        myTree= Tree.get_from_string( treeString, 'newick', annotations_as_nhx=True, extract_comment_metadata=True , suppress_annotations=False)
        return myTree
           
    def checkExe(self, exePath):
        return os.path.isfile(exePath) and os.access(exePath, os.X_OK)
    
    def Where(self, program):
        '''
        input: name of executable
        output: path to executable
        '''
        fpath, fname = os.path.split(program)
        if fpath:
            if self.checkExe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                pathToProgram = os.path.join(path, program)
                if self.checkExe(pathToProgram):
                    return pathToProgram
        return None
    
    def checkNotung(self, program):
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            pathToProgram = os.path.join(path, program)
            if os.path.isfile(pathToProgram):
                return pathToProgram
        return None
    
        
    def WrapFseqboot(self):
        '''
        WrapFseqboot() function is a wrapper over the fseqboot tool.
        input:  inFile     = data file,
                id         = core id (int), 
                seqType    = d (dna); p (protein); r (rna), 
                bootNum    = number of replicates,
                seedNum    = Random number seed between 1 and 32767
                method    = b (Bootstrap) Default, 
                            j (Jackknife)
                            c (Permute species for each character)
                            o (Permute character order)
                            s (Permute within species)
                            r (Rewrite data)),
     
        output: bootstrap alignments to bootstrap_(id).out
        '''
        
        # check if path exists
        cmd= self.Where('fseqboot')
        
        if cmd != None:
            
            bootSeqOutFile= os.path.join(self.work_dir, self.outFile)
            print "bootSeqOutFile: %s"  % bootSeqOutFile
            null = open("/dev/null")
            print "--> Bootstrap begins..."
            try:
                cline = FSeqBootCommandline(sequence = self.inFile, outfile = bootSeqOutFile, seqtype = self.seqType, test= self.method,  seed= self.seedNum, reps = self.bootNum, filter = True)
                stdout, stderr= cline()
                
            except IOError, e:
                print ("Class: Wrapper, Function: WrapFseqboot(): %s " % e)
            print "--> Bootstrap done..."
        else:
            print ("Class: Wrapper, Function: WrapFseqboot(): Error: Path to fseqboot is not set ") 
            sys.exit()

    def WrapFastTreeRooted(self, sTree):
        '''
        WrapFastTreeRooted() function is a wrapper over the fasttree tool.
        
        input:  inputFile: sequences in seqboot output format (Phylip), 
                seqType  : d (dna); p (protein)
                n        :    Number of datasets, 
        
        output: Bifucated and rooted tree in newkick format to file
        '''
        cmd= self.Where('FastTree')
        infile= os.path.join(self.work_dir, self.outFile)
        if cmd != None:
            print "--> Generating newick trees using fasttree begins..."
            try:
                fastreeCmdLine = FastTreeCommandline(cmd=cmd, input=infile, n=self.n, boot=100, nosupport=True, nome= True, seed= self.seedNum, quiet= True)
                child = subprocess.Popen(str(fastreeCmdLine),stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=(sys.platform!="win32"))
                j=1
                with open(os.path.join(self.work_dir, str(self.coreId)+'batch.unrooted.ini'), 'w') as bwf:
                    bwf.write(os.path.abspath(sTree) + '\n')
                    for i in child.stdout:
                        outNewickFile= os.path.join(self.work_dir, str(j)+ "_ftree_coreId_" + str(self.coreId) +".newick")
                        bwf.write(outNewickFile+ '\n')
                        ftreeOutFile= outNewickFile
                        with open(ftreeOutFile, 'w') as wf:
                            
                            tree = Tree.get_from_string(i, 'newick')
                            tree.resolve_polytomies(update_splits=False)
                            st=tree.as_string('newick')
                            #remove the extra strings
                            st='\"'+str(st)+'\"'
                            st=st.translate(None,'\'')
                            st=st.translate(None,'\"')
                            
                            wf.write(st)
                            
                        j += 1

            except IOError, e:
                print ("Class: Wrapper, Function: WrapFastTreeRooted(): %s " % e)
            print "--> Tree generation done..."
        else:
            print ("Class: Wrapper, Function: WrapFastTreeRooted(): Error: Path to 'fasttree' is not set. ") 
            sys.exit()
        
    def reRootGeneTreeUsingNotung(self):
        """
        Root the unrooted output trees from fasttree using Notung
        """
        cmd= self.checkNotung('Notung-2.6.jar')
        ubatchFile= os.path.join(os.path.abspath(self.work_dir),str(self.coreId)+'batch.unrooted.ini')
        if cmd != None:
            print "--> Reconciling gene trees with species tree..."
            try:
                subprocess.call(str("java -jar "+cmd+" -b "+ ubatchFile +" -absfilenames --root --speciestag prefix --edgeweights name --silent --nolosses --treeoutput newick --outputdir "+ os.path.abspath(self.work_dir)), shell=(sys.platform!="win32"), stdout=True, stderr=True)   
            except IOError, e:
                print ("Class: Wrapper, Function: reRootGeneTreeUsingNotung(): %s " % e)
            print "--> Reconciliation done..."
        else:
            print ("Class: Wrapper, Function: reRootGeneTreeUsingNotung(): Error: Path to 'Notung' is not set. ") 
            sys.exit()
    
            
    def WrapMuscle(self, inSeqFile):
        '''
        WrapMuscle() function is a wrapper over the muscle tool.
        input:  inputFile: sequences in fasta format  
        output: MSA to file
        '''
        cmd= self.Where('muscle')
        
        if cmd != None:
            print "--> Generating Perturbed MSA using Muscle begins..."
            for i in xrange(1, self.n+1):
                inputTree=os.path.join(self.work_dir, str(i)+ "_ftree_coreId_" + str(self.coreId) +".newick.rooting.0")
                outFile= os.path.join(self.work_dir, "perturbed_MSA_"+str(i)+".algin")
                try:
                    FNULL = open(os.devnull, 'w')
                    subprocess.call(str(cmd+" -in "+ inSeqFile +" -out "+ outFile + " -usetree_nowarn "+ inputTree), shell=(sys.platform!="win32"), stdout=FNULL, stderr=True)
                    
                except IOError, e:
                    print ("Class: Wrapper, Function: WrapMuscle(): %s " % e)
            print "--> Perturbed MSA generation done..."
        else:
            print ("Class: Wrapper, Function: WrapMuscle(): Error: Path to 'muscle' is not set. ") 
            sys.exit()
            
    def WrapFastTreeUnRooted(self, myFastest, myGtr, myWag, myGama, sTree):
        '''
        WrapFastTreeUnRooted() function is a wrapper over the fasttree tool.
        
        input: None
        
        output: unrooted tree in newkick format to file
        '''
        cmd= self.Where('FastTree')
        infile= os.path.join(self.work_dir, self.outFile)
        if self.seqType == 'p' and myWag == None:
            myWag= True
        elif self.seqType == 'd' and myGtr == None:
            myGtr= True
            
        if cmd != None:
            print "--> Generating newick trees using fasttree begins..."
            try:
                with open(os.path.join(self.work_dir, str(self.coreId)+'batch.unrooted'), 'w') as bwf:
                    with open(os.path.join(self.work_dir, str(self.coreId)+'batch.rooted'), 'w') as rbwf, open(os.path.join(self.work_dir, str(self.coreId)+'batch.resovled'), 'w') as pbwf, open(os.path.join(self.work_dir, str(self.coreId)+'batch.rearranged'), 'w') as rwf:
                        rbwf.write(os.path.abspath(sTree) + '\n')
                        pbwf.write(os.path.abspath(sTree) + '\n')
                        rwf.write(os.path.abspath(sTree) + '\n')
                        bwf.write(os.path.abspath(sTree) + '\n')
                            
                        for i in xrange(1, self.n+1):
                            infile = os.path.join(self.work_dir, "perturbed_MSA_"+str(i)+".algin")
                            outFile= os.path.join(self.work_dir, str(i)+ "_ftree_coreId_"+str(self.coreId)+".newick")
                            fastreeCmdLine = FastTreeCommandline(cmd= cmd, input=infile, quote= True, wag=myWag, boot=100, nosupport=True, nome= True, fastest=myFastest, gtr=myGtr, gamma=myGama, seed= self.seedNum, quiet= True)
                            child = subprocess.Popen(str(fastreeCmdLine),stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=(sys.platform!="win32"))
                            
                            pbwf.write(str(os.path.abspath(str(outFile+".rooting.0.resolve.0")))+ '\n')
                            rbwf.write(str(os.path.abspath(str(outFile+".rooting.0")))+ '\n')
                            rwf.write(str(os.path.abspath(str(outFile+".rooting.0.resolve.0.reconciled")))+ '\n')
                            bwf.write(str(os.path.abspath(outFile))+ '\n')
                            for j in child.stdout:
                                j='\"'+str(j)+'\"'
                                j=j.translate(None,'\'')
                                j=j.translate(None,'\"')
        
                                with open(outFile, 'w') as wf:
                                    wf.write(j)     
            except IOError, e:
                print ("Class: Wrapper, Function: WrapFastTreeUnRooted(): %s " % e)
            print "--> Tree generation done..."
        else:
            print ("Class: Wrapper, Function: WrapFastTreeUnRooted(): Error: Path to 'fasttree' is not set. ") 
            sys.exit()
            
    def RootAndResolvePolytomUsingNotung(self, sTree):
        '''
        This function is a wrapper over the Notung 2.6 tool. It will root an unrooted tree, resovle the polytomies and reconcile with species tree.
        input: gene tree, species tree
        output: reconciled gene tree 
        '''
        cmd= self.checkNotung('Notung-2.6.jar')
        ubatchFile= os.path.join(os.path.abspath(self.work_dir),str(self.coreId)+'batch.unrooted')
        resovledbatchFile= os.path.join(os.path.abspath(self.work_dir),str(self.coreId)+'batch.resovled')
        rbatchFile= os.path.join(os.path.abspath(self.work_dir),str(self.coreId)+'batch.rooted')
        rearrengedBatchFile= os.path.join(os.path.abspath(self.work_dir),str(self.coreId)+'batch.rearranged')
        if cmd != None:
            print "--> Reconciling gene trees with species tree..."
            try:
                
                subprocess.call(str("java -jar "+cmd+" -b "+ ubatchFile +" -absfilenames --root --speciestag prefix --edgeweights name --silent --nolosses --treeoutput newick --outputdir "+ os.path.abspath(self.work_dir)), shell=(sys.platform!="win32"), stdout=True, stderr=True)   
                subprocess.call(str("java -jar "+cmd+" -b "+ rbatchFile +" -absfilenames --resolve --speciestag prefix --edgeweights name --silent --treeoutput newick --outputdir "+ os.path.abspath(self.work_dir)), shell=(sys.platform!="win32"), stdout=True, stderr=True)   
                subprocess.call(str("java -jar "+cmd+" -b "+ resovledbatchFile +" -absfilenames --reconcile --speciestag prefix --edgeweights name --silent --maxtrees 1 --nolosses --outputdir "+ os.path.abspath(self.work_dir)), shell=(sys.platform!="win32"), stdout=True, stderr=True)
                subprocess.call(str("java -jar "+cmd+" -b "+ rearrengedBatchFile +" -absfilenames --rearrange --threshold 50% --speciestag prefix --edgeweights name --silent --maxtrees 1 --nolosses --outputdir "+ os.path.abspath(self.work_dir)), shell=(sys.platform!="win32"), stdout=True, stderr=True)
                
            except IOError, e:
                print ("Class: Wrapper, Function: WrapNotung(): %s " % e)
                
            print "--> Reconciliation done..."
        else:
            print ("Class: Wrapper, Function: WrapNotung(): Error: Path to 'Notung' is not set. ") 
            sys.exit()
    
    
    def bSRapidNJTree(self, sTree):
        '''
        bSRapidNJTree() function is a wrapper over the rapidnj tool.
        
        input:  inputFile: MSA in Stockholm,  
        output: unrooted tree in newkick format
        '''
        cmd= self.Where('rapidnj')
        infile= os.path.join(self.work_dir, self.outFile)
        if cmd != None:
            print "--> Generating newick trees using rapidnj begins..."
            try:
                subprocess.call(str(cmd+" "+ self.inFile +" -i sth -o t -b "+str(self.bootNum)+" -t p >"+ os.path.abspath(self.work_dir)+"/bs.trees"), shell=(sys.platform!="win32"), stdout=True, stderr=True)   
                
                with open(os.path.join(self.work_dir, str(self.coreId)+'batch.rooted'), 'w') as rbwf, open(os.path.join(self.work_dir, str(self.coreId)+'batch.resovled'), 'w') as pbwf, open(os.path.join(self.work_dir, str(self.coreId)+'batch.unrooted'), 'w') as bwf, open(os.path.join(self.work_dir, str(self.coreId)+'batch.rearranged'), 'w') as rwf:
                    rbwf.write(os.path.abspath(sTree) + '\n')
                    pbwf.write(os.path.abspath(sTree) + '\n')
                    rwf.write(os.path.abspath(sTree) + '\n')
                    bwf.write(os.path.abspath(sTree) + '\n')
                    j=1
                    with open(os.path.abspath(self.work_dir)+"/bs.trees", 'r') as rtf:
                        for line in rtf:
                            line= line.strip()
                            line=line.replace("'", "")
                            if not line:
                                pass
                            else:
                                outNewickFile= os.path.join(self.work_dir, str(j)+ "_ftree_coreId_" + str(self.coreId) +".newick")
                                with open(outNewickFile, 'w') as rt:
                                    rt.write(line)
                                rbwf.write(str(os.path.abspath(str(outNewickFile+".rooting.0")))+ '\n')
                                pbwf.write(str(os.path.abspath(str(outNewickFile+".rooting.0.resolve.0")))+ '\n')
                                rwf.write(str(os.path.abspath(str(outNewickFile+".rooting.0.resolve.0.reconciled")))+ '\n')
                                bwf.write(str(os.path.abspath(outNewickFile))+ '\n')
                                j += 1
            except IOError, e:
                print ("Class: Wrapper, Function: bSRapidNJTree(): %s " % e)
            print "--> Tree generation done..."
        else:
            print ("Class: Wrapper, Function: bSRapidNJTree(): Error: Path to 'rapidnj' is not set. ") 
            sys.exit()

    def onlyBSFastTreeRooted(self, myFastest, myGtr, myWag, myGama, sTree):
        '''
        WrapFastTreeRooted() function is a wrapper over the fasttree tool.
        
        input:  inputFile: sequences in seqboot output format (Phylip), 
                seqType  : d (dna); p (protein)
                n        :    Number of datasets, 
        
        output: Bifucated and rooted tree in newkick format to file
        '''
        cmd= self.Where('FastTree')
        infile= os.path.join(self.work_dir, self.outFile)
        if cmd != None:
            print "--> Generating newick trees using fasttree begins..."
            try:
                #fastreeCmdLine = FastTreeCommandline(cmd=cmd, input=infile, n=self.n,quote= True, wag=myWag, noml= True, nome= True, fastest=myFastest, gtr=myGtr, gamma=myGama, seed= self.seedNum, quiet= True)
                fastreeCmdLine = FastTreeCommandline(cmd=cmd, input=infile, n=self.n,quote= True, wag=myWag , boot=100,nome= True, fastest=myFastest, gtr=myGtr, gamma=myGama, seed= self.seedNum, quiet= True)
                child = subprocess.Popen(str(fastreeCmdLine),stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=(sys.platform!="win32"))
                j=1
                with open(os.path.join(self.work_dir, str(self.coreId)+'batch.rooted'), 'w') as rbwf, open(os.path.join(self.work_dir, str(self.coreId)+'batch.resovled'), 'w') as pbwf, open(os.path.join(self.work_dir, str(self.coreId)+'batch.unrooted'), 'w') as bwf, open(os.path.join(self.work_dir, str(self.coreId)+'batch.rearranged'), 'w') as rwf:
                    rbwf.write(os.path.abspath(sTree) + '\n')
                    pbwf.write(os.path.abspath(sTree) + '\n')
                    rwf.write(os.path.abspath(sTree) + '\n')
                    bwf.write(os.path.abspath(sTree) + '\n')
                    for i in child.stdout:
                        #print "tree %d" % j
                        outNewickFile= os.path.join(self.work_dir, str(j)+ "_ftree_coreId_" + str(self.coreId) +".newick")
                        ftreeOutFile= outNewickFile
                        with open(ftreeOutFile, 'w') as wf:
                            i='\"'+str(i)+'\"'
                            i=i.translate(None,'\'')
                            i=i.translate(None,'\"')
                            wf.write(i)
#                            tree = Tree.get_from_string(i, 'newick')
#                            #tree.resolve_polytomies(update_splits=False)
#                            st=tree.as_string('newick')
#                            wf.write(st)
                            pbwf.write(str(os.path.abspath(str(outNewickFile+".rooting.0.resolve.0")))+ '\n')
                            rbwf.write(str(os.path.abspath(str(outNewickFile+".rooting.0")))+ '\n')
                            rwf.write(str(os.path.abspath(str(outNewickFile+".rooting.0.resolve.0.reconciled")))+ '\n')
                            bwf.write(str(os.path.abspath(outNewickFile))+ '\n')
                        j += 1

            except IOError, e:
                print ("Class: Wrapper, Function: WrapFastTreeRooted(): %s " % e)
            print "--> Tree generation done..."
        else:
            print ("Class: Wrapper, Function: WrapFastTreeRooted(): Error: Path to 'fasttree' is not set. ") 
            sys.exit()
            
    def onlyBSReconcileUsingNotung(self, sTree):
        '''
        This function is a wrapper over the Notung 2.6 tool. It will root an unrooted tree, resovle the polytomies and reconcile with species tree.
        input: gene tree, species tree
        output: reconciled gene tree 
        '''
        cmd= self.checkNotung('Notung-2.6.jar')
        ubatchFile= os.path.join(os.path.abspath(self.work_dir),str(self.coreId)+'batch.unrooted')
        resovledbatchFile= os.path.join(os.path.abspath(self.work_dir),str(self.coreId)+'batch.resovled')
        rbatchFile= os.path.join(os.path.abspath(self.work_dir),str(self.coreId)+'batch.rooted')
        rearrengedBatchFile= os.path.join(os.path.abspath(self.work_dir),str(self.coreId)+'batch.rearranged')
        if cmd != None:
            print "--> Reconciling gene trees with species tree..."

            try:
                subprocess.call(str("java -jar "+cmd+" -b "+ ubatchFile +" -absfilenames --root --speciestag prefix --edgeweights name --silent --nolosses --treeoutput newick --outputdir "+ os.path.abspath(self.work_dir)), shell=(sys.platform!="win32"), stdout=True, stderr=True)   
                subprocess.call(str("java -jar "+cmd+" -b "+ rbatchFile +" -absfilenames --resolve --speciestag prefix --edgeweights name --silent --treeoutput newick --outputdir "+ os.path.abspath(self.work_dir)), shell=(sys.platform!="win32"), stdout=True, stderr=True)   
                subprocess.call(str("java -jar "+cmd+" -b "+ resovledbatchFile +" -absfilenames --reconcile --speciestag prefix --edgeweights name --silent --maxtrees 1 --nolosses --outputdir "+ os.path.abspath(self.work_dir)), shell=(sys.platform!="win32"), stdout=True, stderr=True)
                subprocess.call(str("java -jar "+cmd+" -b "+ rearrengedBatchFile +" -absfilenames --rearrange --threshold 50% --speciestag prefix --edgeweights name --silent --maxtrees 1 --nolosses --outputdir "+ os.path.abspath(self.work_dir)), shell=(sys.platform!="win32"), stdout=True, stderr=True)
                
            except IOError, e:
                print ("Class: Wrapper, Function: WrapNotung(): %s " % e)
                
            print "--> Reconciliation done..."
        else:
            print ("Class: Wrapper, Function: WrapNotung(): Error: Path to 'Notung' is not set. ") 
            sys.exit()
            
    def reconcileGtreeUsingNotung(self, sTree, gTree):
        '''
        This function is a wrapper over the Notung 2.6 tool. It will root an unrooted tree, resovle the polytomies and reconcile with species tree.
        input: gene tree, species tree
        output: reconciled gene tree 
        '''
        cmd= self.checkNotung('Notung-2.6.jar')
        #rOutFile= os.path.join(os.path.abspath(self.work_dir),str(self.coreId)+'batch.rooted')
        if cmd != None:
            print "--> Reconciling gene trees with species tree..."  
            try:
                #subprocess.call(str("java -jar "+cmd+" -b "+ batchFile +" -absfilenames --root --speciestag prefix  --nolosses --treeoutput newick --outputdir "+ os.path.abspath(self.work_dir)), shell=(sys.platform!="win32"), stdout=True, stderr=True)   
                #gTree= gTree+".rooting.0"
                subprocess.call(str("java -jar "+cmd+" -g "+ gTree +" -s "+ sTree + " --reconcile --speciestag prefix --edgeweights name --silent --maxtrees 1 --nolosses --outputdir "+ self.work_dir), shell=(sys.platform!="win32"), stdout=True, stderr=True)
                
            except IOError, e:
                print ("Class: Wrapper, Function: WrapNotung(): %s " % e)            
            print "--> Reconciliation done..."
        else:
            print ("Class: Wrapper, Function: WrapNotung(): Error: Path to 'Notung' is not set. ") 
            sys.exit()
            
    
    def setSpeciesLeafLabels(self, streeString):
        '''
            input: path to species tree
        '''
        sTree= self.readTreeFromString(streeString)

        self.setSpeceistree(sTree)
        self.leafLabelStree, self.internalExRootSpeceLabels = self.allLeafLabels(sTree)
    
    def cutSubStringFromNotungSTree(self, myString, nbits):
        '''
        Remove the first nbits from the string and the last ']'
            
        '''    
        myString=myString[nbits:]
        myString=myString[:-1]+';'
        
        return myString   
        
    def readRecGtreeAndCut(self):
        '''
            This function takes all the reconcilizations and cut at speciation nodes
        '''
        orthoFamilies       = {} # key= root label of subtree; value: subtree
        paraFamailies       = {}
        speciesOrthoCuts    = {}
        speciesParaCuts     = {}
        
        #with open(os.path.join(os.path.abspath(self.work_dir),str(self.coreId)+'batch.rearranged') , 'r') as rf:
        with open(os.path.join(os.path.abspath(self.work_dir),str(self.coreId)+'batch.resovled') , 'r') as rf:
            readFlag=0
            for line in rf:
                if readFlag != 0:
                    
                    line=line.rstrip()
                    #linesTree=0
                    # check if file exists
                    #file=line +'.rearrange'
                    file=line +'.reconciled'
                    if os.path.isfile(file) == True:
                        pass
                    else:
                        file=file+'.0'
                    with open(file, 'r') as f:
                        readStreeFlag= True
                        for t in f:
                            
                            t=t.rstrip()
                            recGTree= self.readTreeFromString(t)
                            if readStreeFlag == False:
                                break
                            else:
                                t= f.next()
                                t=t.rstrip()
                                t= self.cutSubStringFromNotungSTree(t, 22)
                                #print t
                                self.setSpeciesLeafLabels(t)
                                readStreeFlag = False
                                break
                                
                        self.cutAtSpeciationVertices( orthoFamilies, paraFamailies, recGTree, speciesOrthoCuts, speciesParaCuts)

                    readFlag +=1
                else:
                    
                    readFlag = 1 
            try:
                
                orthoHandle = csv.writer(open(self.resultsFile+"."+str(self.coreId)+".sallfamilies",  'wb'))
                
                for key, value in orthoFamilies.items():
                    orthoHandle.writerow([key, value, speciesOrthoCuts[key]  ])

                paraHandle = csv.writer(open(self.resultsFile+"."+str(self.coreId)+".dallfamilies",  'wb'))
                for key, value in paraFamailies.items():
                    paraHandle.writerow([key, value, speciesParaCuts[key] ])

            except IOError, e:
                print ("Class: Wrapper, Function: readRecGtreeAndCut(): %s " % e)

        
    def cutAtSpeciationVertices(self, orthoFamilies, paraFamailies , recGTree, speciesOrthoCuts, speciesParaCuts):
        '''
            This function cuts the reconciled gene tree at the speciation vertices.
            input:     Reconciled Gene Tree, path to the file containing species tree
            output:    families. soft orthologes and soft paraloges
        '''
        
        for i in recGTree.internal_nodes():
            leafLabels      =[]

            d= i._annotations.values_as_dict()['D']
            s= i._annotations.values_as_dict()['S']
            if  d == 'N':
                # orthologous families
                t= Tree(i)
                leafLabels= self.onlyLeafLabels(t)
                family=""
                family= str(sorted(leafLabels))
                
                if s not in speciesOrthoCuts.keys():
                    speciesOrthoCuts[family]=s
                if family not in orthoFamilies.keys():
                    orthoFamilies[family]=1
                else:
                    orthoFamilies[family] +=1
            elif s in self.leafLabelStree:
                # paralogous families
                t= Tree(i)
                leafLabels= self.onlyLeafLabels(t)
                family=""
                family= str(sorted(leafLabels))
                if s not in speciesParaCuts.keys():
                    speciesParaCuts[family]=s
                if family not in paraFamailies.keys():
                    paraFamailies[family] =1
                else:
                    paraFamailies[family] +=1

            elif s in self.internalExRootSpeceLabels:
                # paralogous families
                t= Tree(i)
                leafLabels= self.onlyLeafLabels(t)
                family=""
                family= str(sorted(leafLabels))
                if s not in speciesParaCuts.keys():
                    speciesParaCuts[family]=s
                if family not in paraFamailies.keys():
                    paraFamailies[family] = 1
                else:
                    paraFamailies[family] +=1
                
                

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
