#!/usr/bin/env python
# encoding: utf-8
"""
TriggerStimulusGenerator.py

Created by  on 2006-08-03.
Copyright (c) 2006 Barry Wark. All rights reserved.
"""

from Foundation import *
from acqui.AcqUI import *
from numpy import *

from math import ceil

import unittest
import objc


class TriggerStimulusGenerator (AUIStimulusPyBase):
    def init(self):
        self = super(TriggerStimulusGenerator, self).init()
        
        if(self != None):
            self.triggerBits = []
            self.triggerOnSamples = []
            self.dtypeStr = dtypeForDigitalOutput().str #array((), dtype=int_).dtype.str
        
        return self
    
        
    def stimulusData(self):
        return self.stimulusDataForParamsDictionary_(self.paramsDictionary())
        
    
        
    def stimulusDataForParamsDictionary_version_(self, p, v):
        triggerBits = p['triggerBits']
        triggerOnSamples = p['triggerOnSamples']
        dataType = dtype(p['dtypeStr'])
        
        samples = ceil(self.length * self.sampleRate)
        
        arr = zeros(samples, dtype=dataType)
        triggerSample = 0
        for b in self.triggerBits:
            t = 1 << b
            triggerSample += t
        
        arr[self.triggerOnSamples] = triggerSample
        
        return arr
    
    
    def paramsDictionary(self):
        d = dict(
            triggerOnSamples = self.triggerOnSamples,
            triggerBits = NSArray.arrayWithArray_(self.triggerBits),
            dtypeStr = self.dtypeStr
        )
        
        return NSDictionary.dictionaryWithDictionary_(d)
    
        
    
    
    def version(self):
        return 1
    
    
    @classmethod
    def stimulusID(self):
        return 'edu.washington.bwark.acqui.stim.trigger'        
    

objc.removeAutoreleasePool()

class TriggerStimulusGeneratorTests(unittest.TestCase):
    def setUp(self):
        pass
    
    def testTrigger(self):
        stim = TriggerStimulusGenerator.alloc().init()
        
        for i in xrange(10):
            stim.sampleRate = 10000
            stim.length = 1
            stim.triggerAmp = random.rand()
            tSamples = floor(random.rand(100)*10000).astype(int_)
        
            s = stim.stimulusData()
            d = stim.paramsDictionary()
            t = zeros(ceil(stim.length * stim.sampleRate))
            t[tSamples] = t.triggerAmp
        
            self.assertTrue(all(s == t))
            self.assertTrue(all(s[tSamples] == stim.triggerAmp))
            self.assertTrue(all(s == stim.stimulusDataForParamsDictionary_(d)))
    


if __name__ == '__main__':
    unittest.main()
