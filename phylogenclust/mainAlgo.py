'''
Created on Sep 19, 2013

@author: Mehmood Alam Khan Malagori
@email:   malagori@kth.se
'''
__version__= "0.1dev"
__license__ = "GPLv3"



import random as rn
from Bio import SeqIO
import os

from wrappers import Wrapper

class MainAlgo(object):
    '''
    This class contain the main work flow of the algorithm.
    '''


    def __init__(self, outDir, resultsFile, dataFile, fileFormat, seqType, bootNum, method, interLeaved, fastest, gtr, wag, gama, rSeed, sTree, noMSA, rapidNJFlag, myRank):
        '''
        Constructor
        '''
        self.randomSeed         = rSeed
        self.workDir            = outDir
        self.dataFile           = dataFile
        self.seqType            = seqType
        self.bootNum            = bootNum
        self.resultsFile        = resultsFile
        self.samplingMethod     = method
        self.interLeaved        = interLeaved
        self.fileFormat         = fileFormat
        self.fastaFile          =''
        self.phylipFile         =''
        self.stockholmFile      =''
        self.coreId             = myRank
        self.fastest            = fastest
        self.gtr                = gtr
        self.wag                = wag
        self.gama               = gama
        self.sTree              = sTree
        self.noMSA              = noMSA
        self.rapidNJFlag        = rapidNJFlag
 
        
    def touch(self,path):
        '''
            create empty file
        '''
        with open(path, 'a'):
            os.utime(path, None)
        
    def convertFastaToPhylip(self, infile):
        '''
            This function will convert fasta file to phylip format
            input: Path to file containing fasta sequences
            output: Path to file containing phylip format sequences
        '''
        
        try:
            self.phylipFile= self.workDir+"/inputMSA.phylip"
            
            self.touch(self.phylipFile)
                
            count = SeqIO.convert(infile, "fasta", self.phylipFile, "phylip")
            
        except IOError, e:

            print ("Class: MainAlgo, Function: convertFastaToPhylip(): %s " % e)
 
    def convertPhylipToFasta(self, infile):
        '''
            This function will convert phylip file to fasta format
            input: Path to file containing phylip sequences
            output: Path to file containing fasta format sequences  
        '''
        try:
            self.fastaFile= self.workDir+"/inputMSA.fasta"
            self.touch(self.fastaFile)
            count = SeqIO.convert(infile, "phylip", self.fastaFile, "fasta")

        except IOError, e:
            print ("Class: MainAlgo, Function: convertPhylipToFasta(): %s " % e)
    
    def convertPhylipToStockholm(self, infile):
        '''
            This function will convert phylip file to stockholm format
            input: Path to file containing phylip sequences
            output: Path to file containing stockholm format sequences  
        '''
        try:
            self.stockholmFile= self.workDir+"/inputMSA.stockholm"
            self.touch(self.stockholmFile)
            count = SeqIO.convert(infile, "phylip", self.stockholmFile, "stockholm")

        except IOError, e:
            print ("Class: MainAlgo, Function: convertPhylipToStockholm(): %s " % e)
            
    def convertFastaToStockholm(self, infile):
        '''
            This function will convert fasta file to stockholm format
            input: Path to file containing fasta sequences
            output: Path to file containing stockholm format sequences  
        '''
        try:
            self.stockholmFile= self.workDir+"/inputMSA.stockholm"
            self.touch(self.stockholmFile)
            count = SeqIO.convert(infile, "fasta", self.stockholmFile, "stockholm")

        except IOError, e:
            print ("Class: MainAlgo, Function: convertFastaToStockholm(): %s " % e)
    
    def runAlgo(self):
        '''
        This function controls the main work flow of the GFI algorithm.
        '''
        # convert the input file to required file format
        if self.fileFormat == 'F':
            self.fastaFile= self.dataFile
            self.convertFastaToPhylip(self.dataFile)
        elif self.fileFormat == 'P':
            self.phylipFile= self.dataFile
            self.convertPhylipToFasta(self.dataFile)
        
        # if rapidnj is used as tree reconstruction program, then
        if self.rapidNJFlag == True:
            
            # generate stockholm format
            self.convertPhylipToStockholm(self.phylipFile)
            # create Wrapper class object
            wrapperObj= Wrapper(self.workDir, self.resultsFile, self.stockholmFile, self.coreId, self.seqType, self.randomSeed, self.bootNum, self.samplingMethod, self.interLeaved)
            # generate bootstrap trees
            wrapperObj.bSRapidNJTree(self.sTree)
            # reconcile bs trees with species tree
            wrapperObj.onlyBSReconcileUsingNotung(self.sTree)
        else:
            
            # instantiate Wrapper class object
            wrapperObj= Wrapper(self.workDir, self.resultsFile, self.phylipFile, self.coreId, self.seqType, self.randomSeed, self.bootNum, self.samplingMethod, self.interLeaved)
            
            # call WrapFseqboot Fucntion
            wrapperObj.WrapFseqboot()
            
            if self.noMSA == True: # disabling the perturbation step
                wrapperObj.onlyBSFastTreeRooted(self.fastest, self.gtr, self.wag, self.gama, self.sTree)
                wrapperObj.onlyBSReconcileUsingNotung(self.sTree)
            else:
                # call WrapFastTree Function
                wrapperObj.WrapFastTreeRooted(self.sTree)
                # call the reRoot function
                #wrapperObj.reRootGeneTreeUsingNotung()
                # call WrapMuscle Function
                wrapperObj.WrapMuscle(self.fastaFile)
                # call WrapFastTree Function
                wrapperObj.WrapFastTreeUnRooted(self.fastest, self.gtr, self.wag, self.gama, self.sTree)
                # call RootAndResolvePolytomUsingNotung Function
                wrapperObj.RootAndResolvePolytomUsingNotung(self.sTree)
            
        # call readRecGtreeAndCut function
        wrapperObj.readRecGtreeAndCut()
            
            #with open(self.resultsFile+ ".speciestree", 'w') as wf:
            #    wf.write(wrapperObj.getSpeciesTree())
        # this species tree is the species tree from the last bootstrap gene tree reconciliation phase.
        # in order to get the whole tree you can print the original (big) species tree from the notung to get all the internal vertices of species tree.
        wrapperObj.writeSpeciesTree()
        
        
        
        
        