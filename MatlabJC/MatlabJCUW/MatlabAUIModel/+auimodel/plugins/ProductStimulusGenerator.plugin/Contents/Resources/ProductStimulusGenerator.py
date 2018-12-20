#!/usr/bin/env python
# encoding: utf-8
"""
ProductStimulusGenerator.py

Created by  on 2007-03-05.
Copyright (c) 2007 Barry Wark. All rights reserved.
"""

#from Foundation import *
from acqui.AcqUI import AUICompoundStimulusGenerator

import numpy as N

import objc

class ProductStimulusGenerator (AUICompoundStimulusGenerator):
    
    def init(self):
        self = super(ProductStimulusGenerator, self).init()
        
        if(self != None):
            self.combFn = lambda (s): N.prod(N.vstack(s), axis=0)
        
        return self
    
    
    def version(self):
        return 1
     
    @classmethod
    def stimulusID(self):
        return 'edu.washington.bwark.acqui.stim.product'        
    

objc.removeAutoreleasePool()