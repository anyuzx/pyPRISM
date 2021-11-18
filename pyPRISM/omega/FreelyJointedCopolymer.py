#!python
from pyPRISM.omega.Omega import Omega
import numpy as np


class FreelyJointedCopolymer(Omega):
    r'''Freely jointed chain intra-molecular correlation function for copolymer
    
    **Description**
        
        the freely-jointed chian is an ideal polymer chain model that assumes
        a constant bond length :math:`l` and no correlations between the
        directions of different bond vectors. If the chain is a copolymer
        with two types of monomer. The correlation function can be computed
        between different pairs.
    '''
    
    def __init__(self, sequence, l, pair1, pair2):
        r'''Constructor
        
        Arguments
        ---------
        sequence: string
            a string containing only two types of string
        
        l: float
            bond length
            
        pair: tuple
            tuple of pair
        '''
        self.sequence = np.array(list(sequence))
        self.length = self.N = len(sequence)
        self.l = l
        self.value = None
        self.pair1 = pair1
        self.pair2 = pair2
        
        assert self.pair1 in np.unique(self.sequence) and self.pair2 in np.unique(self.sequence), "Pairs specified are not in the sequence"
        
    
    def __repr__(self):
        return '<Omega: FreelyJointedCopolymer>'
    
    def calculate(self, k):
        '''Return value of :math:`\hat{\omega}` at supplied :math:`k`
        
        Arguments
        ---------
        k: np.ndarray
            array of wavenumber values to calculate :math:`\omega` at
        
        '''        
        if self.pair1 == self.pair2:
            index = np.where(self.sequence == self.pair1)[0]
            index_num = index.shape[0]
            
            self.value = np.sum(np.power((np.sin(k * self.l) / (k * self.l))[:, np.newaxis, np.newaxis], np.abs(index[:, np.newaxis] - index)[np.newaxis,:]), axis=(1,2)) / (index_num)
        else:
            index_1 = np.where(self.sequence == self.pair1)[0]
            index_2 = np.where(self.sequence == self.pair2)[0]
            
            index_1_num = index_1.shape[0]
            index_2_num = index_2.shape[0]
            self.value = np.sum(np.power((np.sin(k * self.l) / (k * self.l))[:, np.newaxis, np.newaxis], np.abs(index_1[:, np.newaxis] - index_2)[np.newaxis,:]), axis=(1,2)) / (index_1_num + index_2_num)
        
        return self.value
