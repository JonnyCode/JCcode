#!/usr/bin/env python
# encoding: utf-8
"""
SumStimulusGenerator.py

Created by  on 2006-08-03.
Copyright (c) 2006 Barry Wark. All rights reserved.
"""

from Foundation import NSLog
from acqui.AcqUI import AUICompoundStimulusGenerator
import numpy as N

import unittest
import objc


class SumStimulusGenerator (AUICompoundStimulusGenerator):
    
    def init(self):
        self = super(SumStimulusGenerator, self).init()
            
        if(self != None):
            self.combFn = lambda(s): N.sum(N.vstack(s), axis=0)
        
        return self
    
    
    def version(self):
        return 1
    
    @classmethod
    def stimulusID(self):
        return 'edu.washington.bwark.acqui.stim.sum'
    


objc.removeAutoreleasePool()


if __name__ == '__main__':
    unittest.main()
