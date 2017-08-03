import unittest
import numpy as np
import typyPRISM

class PRISM_TestCase(unittest.TestCase):
    def setup(self):
        '''Construct a simple but fully specified PRISM problem'''
        
        sys = typyPRISM.System(['A','B'])
        
        sys.domain = typyPRISM.Domain(dr=0.1,length=512)
        
        sys.density['A'] = 0.5
        sys.density['B'] = 0.5
        
        sys.closure.setUnset(typyPRISM.closure.PercusYevick())
        
        sys.potential.setUnset(typyPRISM.potential.HardSphere(sigma=1.0))
        
        sys.intraMolCorr['A','A'] = typyPRISM.intraMolCorr.SingleSite()
        sys.intraMolCorr['A','B'] = typyPRISM.intraMolCorr.NoIntra()
        sys.intraMolCorr['B','B'] = typyPRISM.intraMolCorr.Gaussian(sigma=1.0,length=10)
        
        PRISM = sys.createPRISM()
        
        return PRISM
    def test_solve(self):
        '''Can we solve the PRISM equations?'''
        PRISM = self.setup()
        result = PRISM.solve(disp=False)
        self.assertIsNot(result,None)
        
        
        