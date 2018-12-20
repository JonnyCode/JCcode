#!/usr/bin/env python
# encoding: utf-8
"""
StepStimulusGenerator.py

Created by  on 2006-08-03.
Copyright (c) 2006 Barry Wark. All rights reserved.
"""

from Foundation import *
from acqui.AcqUI import *
from numpy import *

from math import ceil

import unittest
import objc


class StepStimulusGenerator (AUIStimulusPyBase):
    def init(self):
        self = super(StepStimulusGenerator, self).init()
        
        if(self != None):
            self.lowTime = 0
            self.highTime = 0
            self.postTime = 0
            self.lowMean = 0
            self.highMean = 0
            self.nReps = 1
            self.sampleRate = 0
        
        return self
    
        
    def stimulusData(self):
        return self.stimulusDataForParamsDictionary_(self.paramsDictionary())
        
    
        
    def stimulusDataForParamsDictionary_version_(self, p, v):
        lowTime = p['lowTime']
        highTime = p['highTime']
        postTime = p.get('postTime', 0) # version<3 d/n have postTime
        lowMean = p['lowMean']
        highMean = p['highMean']
        nReps = p['nReps']
        sampleRate = p.get('sampleRate', 10000) #old versions didn't save sampleRate
        
        lowSamples = ceil(sampleRate * lowTime)
        highSamples = ceil(sampleRate * highTime)
        postSamples = ceil(sampleRate * postTime)
        
        arr = r_[zeros(lowSamples, dtype=float_)+lowMean, zeros(highSamples, dtype=float_)+highMean, 
                zeros(postSamples, dtype=float_)+lowMean]
        arr = tile(arr, (1, nReps))
        
        return arr.astype(dtypeForAnalogOutput())
    
    
    def paramsDictionary(self):
        d = dict(
            lowTime = self.lowTime,
            highTime = self.highTime,
            postTime = self.postTime,
            lowMean = self.lowMean,
            highMean = self.highMean,
            nReps = self.nReps,
            sampleRate = self.sampleRate,
        )
        
        return NSDictionary.dictionaryWithDictionary_(d)
    
        
    
    
    def version(self):
        return 3
    
    
    @classmethod
    def stimulusID(self):
        return 'edu.washington.bwark.acqui.stim.step'        
    

objc.removeAutoreleasePool()

class StepStimulusGeneratorTests(unittest.TestCase):
    def setUp(self):
        pass
    
    def testStepSize(self):
        stim = StepStimulusGenerator.alloc().init()
        
        for i in xrange(10):
            stim.preTime = numpy.random.rand()
            stim.postTime = numpy.random.rand()
            stim.stepTime = numpy.random.rand()
            
            stim.startAmp = numpy.random.rand()
            stim.stepAmp = stim.startAmp + numpy.random.rand()
            stim.sampleRate = 10000 * numpy.random.rand()
            
            s = stim.stimulusData()
            d = stim.paramsDictionary()
            
            self.assertEquals(s.size, (stim.preTime+stim.PostTime+stim.stepTime)/stim.sampleRate)
            self.assertTrue(numpy.all(s == stim.stimulusDataForParamsDictionary(d)))
    
    
    def testStep(self):
        stim = StepStimulusGenerator.alloc().init()
        
        stim.sampleRate = 10000
        stim.preTime = 0.1
        stim.stepTime = 0.1
        stim.postTime = 0.1
        stim.startAmp = 0
        stim.stepAmp = 1
        
        s = stim.stimulusData()
        
        self.assertTrue(numpy.all(s == r_[zeros(1000), zeros(1000)+1, zeros(1000)]))
    


if __name__ == '__main__':
    unittest.main()
