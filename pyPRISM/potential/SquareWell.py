#!python
from __future__ import division,print_function
from pyPRISM.potential.Potential import Potential
import numpy as np

class SquareWell(Potential):
    r'''Square well attractive or repulsive interactions
    
    
    **Mathematical Definition**
    
    .. math::

        U_{\alpha,\beta}(r \geq \sigma_{\alpha,\beta} and r \leq \sigma_{\alpha,\beta} + \alpha_{\alpha,\beta}) = \epsilon_{\alpha,\beta} 

    .. math::

        U_{\alpha,\beta}(r < \sigma_{\alpha,\beta}) = C^{high}

    
    **Variable Definitions**
    
    :math:`\alpha_{\alpha,\beta}`
        Width of square well potential between sites
    :math:`\alpha` and :math:`\beta`.

    :math:`\sigma_{\alpha,\beta}`
        Contact distance of interactions between sites 
	:math:`\alpha` and :math:`\beta`.


    :math:`\epsilon_{\alpha,\beta}`
        Interaction strength between sites 
	:math:`\alpha` and :math:`\beta`.

    :math:`C^{high}`
        High value used to approximate an infinite potential due to overlap
    

    **Description**

        This potential models a square-well form of attraction or repulsion
        between sites with a specified site size and contact distance.

    Example
    -------
    .. code-block:: python

        import pyPRISM
	
        #Define a PRISM system and set the A-B interaction potential
        sys = pyPRISM.System(['A','B'],kT=1.0)
        sys.domain = pyPRISM.Domain(dr=0.1,length=1024)
        sys.potential['A','B'] = pyPRISM.potential.SquareWell(epsilon=1.0,sigma=8.0,alpha=0.5,high_value=10**6)

    .. warning::

        If sigma is specified such that it does not fall on the solution grid
        of the :class:`~pyPRISM.core.Domain.Domain` object specified in
        :class:`~pyPRISM.core.System.System`, then the sigma will effectively
        be rounded. A warning should be emitted during the construction of a
        :class:`~pyPRISM.core.PRISM.PRISM` object if this occurs.
    
    '''
    def __init__(self,epsilon,alpha,sigma=None,high_value=1e6):
        r''' Constructor
        
        Arguments
        ---------
        epsilon: float
            Strength of attraction (positive) or repulsion (negative)
            
        alpha: float
            Range of attraction

        sigma: float, *optional*
            Contact distance. If not specified, sigma will be calculated from 
            the diameters specified in the :class:`~pyPRISM.core.System.System`
            object.
            
        high_value: float, *optional*
            High value used to approximate an infinite potential due to overlap
        
        '''
        self.epsilon = epsilon
        self.sigma = sigma
        self.alpha = alpha
        self.high_value = high_value
        self.funk  = lambda r,sigma: - epsilon
    def __repr__(self):
        return '<Potential: Exponential>'
    
    def calculate(self,r):
        r'''Calculate the value of the potential

        Attributes
        ----------
        r: float np.ndarray
            Array of pair distances at which to calculate potential values
        '''
        assert (self.sigma is not None), 'Sigma must be set before evaluating potential!'

        magnitude = self.funk(r,self.sigma)
        magnitude = np.where(r<=self.sigma, self.high_value, magnitude)
        magnitude = np.where(r>=self.sigma + self.alpha, 0.0, magnitude)
        return magnitude
      
