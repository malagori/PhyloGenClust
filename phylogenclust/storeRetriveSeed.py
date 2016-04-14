'''
Created on Sep 19, 2013

@author: Mehmood Alam Khan Malagori
@email:   malagori@kth.se
'''
__version__= "0.1dev"
__license__ = "GPLv3"


import random
import os
import cPickle as pickle

class RandomSeed(object):
    def __init__(self, infile=None):
        self.state=None
        

    def getSateFromFile(self,infile):
        """ input: infile, File containing the state of random number generator function
            Return: state for random generator function stored in a file. 
        """
        if infile != None:
            if os.path.exists(infile):
                # Restore the previously saved sate
                print '%s, initializing random module' % infile
                with open(infile, 'rb') as f:
                    state = pickle.load(f)
                self.setState(state)
                return self.getState()
            else:
                # Use a well-known start state
                print 'No seed file exist'
    
    def getState(self):
        return self.state
    
    def setState(self,initialState):
        self.state= initialState
        
    def setInitialState(self,seed):
        """
        Args: initialSeed 
        Return: None. but set the state with initialSeed provided
        """
        random.seed(seed)
        self.setState(random.getstate())
        
    def storeSate(self,infile):
        # Save state to a file for next time usage
        # file name should include: dataset name, iteration number, initial seed provided by user
        with open(infile, 'wb') as f:
            pickle.dump(random.getstate(), f)
        