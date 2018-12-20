#!/usr/bin/env python
# encoding: utf-8
"""
SineStimulusGenerator.py

Created by  on 2007-02-17.
Copyright (c) 2007 Barry Wark. All rights reserved.
"""

from __future__ import division

from acqui import AcqUI
import numpy as np

import objc

import unittest

class SineStimulusGenerator(AcqUI.AUIStimulusPyBase):
    def init(self):
        self = super(SineStimulusGenerator, self).init()
        
        if(self != None):
            self.mean = 0
            self.frequency = 0
            self.amp = 0
            self.startTime = 0
            self.phase = 0
        
        return self
    
        
    def stimulusDataForParamsDictionary_version_(self, p, v):
        mean = float(p['mean'])
        phase = float(p['phase'])
        f = float(p['frequency']) #in Hz
        A = float(p['amp'])
        sampleRate = float(p['sampleRate'])
        length = float(p['length'])
        dt = float(p['startTime'])
        
        assert(f <= sampleRate/2) #make sure there won't be aliasing
        
        time = np.linspace(0, length, num=round(sampleRate*length))
        
        arr = (A * np.sin(2*np.pi*f*(time + dt) - phase)) + mean
        
        return arr.astype(AcqUI.dtypeForAnalogOutput())
    
    
    def description(self):
        return "SineStimulusGenerator: " + str(self.paramsDictionary())
    
    
    def paramsDictionary(self):
        d = dict(
            mean = self.mean,
            frequency = self.frequency,
            phase = self.phase,
            amp = self.amp,
            sampleRate = self.sampleRate,
            length = self.length,
            startTime = self.startTime
        )
        
        return AcqUI.cocoaCollection(d)
    
        
    
    
    def version(self):
        return 1
    
    @classmethod
    def stimulusID(self):
        return 'edu.washington.bwark.acqui.stim.sine'        
    


objc.removeAutoreleasePool()


class SineStimulusGeneratorTests(unittest.TestCase):
    def setUp(self):
        self.stim = SineStimulusGenerator.alloc().init()
        
        self.stim.length = 1
        self.stim.sampleRate = 10000
        
        self.stim.mean = 0.0
        self.stim.frequency = 100.0
        self.stim.phase = np.pi
        self.stim.amp = 0.2
    
    
    def testSine(self):
        stim = self.stim.stimulusData()
        t = np.linspace(0,self.stim.length, self.stim.length*self.stim.sampleRate)
        self.assertTrue(np.all(stim == self.stim.amp*np.sin((2*np.pi*t) - phase)+mean))
    
    
    def testThrowsForAliasing(self):
        self.stim.sampleRate = 1
        self.stim.frequency = 100.0
        self.assertThrows(self.stim.stimulusData())
    


if __name__ == '__main__':
    unittest.main()
