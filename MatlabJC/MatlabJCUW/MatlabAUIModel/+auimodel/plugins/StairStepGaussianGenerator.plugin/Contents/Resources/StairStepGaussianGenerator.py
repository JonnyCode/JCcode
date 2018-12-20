#!/usr/bin/env python
# encoding: utf-8
"""
StairStepGaussianGenerator.py

Created by Barry Wark on 2007-07-19.
Copyright (c) 2007 Barry Wark. All rights reserved.
"""

from __future__ import division

from acqui import AcqUI
from Foundation import NSLog

import acqui.noise_source as noise
import numpy as N
import cPickle as pickle

import objc

class StairStepGaussianGenerator(AcqUI.AUIStimulusPyBase):
    
    def init(self):
        self = super(StairStepGaussianGenerator, self).init()
        
        if(self != None):
            self.mean = 0
            self.var = 1
            self.stepRate = 1
            self.length = 1
            self.sampleRate = 1e4
        
        return self
    
    
    
    def stimulusDataForParamsDictionary_version_(self, p, v):
        mean = float(p['mean'])
        var = float(p['var'])
        sampleRate = float(p['sampleRate'])
        length = float(p['length'])
        stepRate = float(p['stepRate'])
        
        if(N.__version__ != p['numpyVersion']):
            NSLog('Warning: A different numpy version (' + p['numpyVersion'] + ') was used to create SSG stimulus. Current version: ' + N.__version__)
            
        self.rSource.set_state(pickle.loads(str(p['randomState'])))
        
        ns = noise.StairStepNoiseSource(noise.GaussianNoiseSource(mean, var, self.rSource), stepRate)
        
        arr = ns.noise(length, sampleRate)
        
        return arr.astype(AcqUI.dtypeForAnalogOutput())
    
    
    def description(self):
        return "StairStepGuassianGenerator: " + str(self.paramsDictionary())
    
    
    def paramsDictionary(self):
        d = dict(
            mean = self.mean,
            var = self.var,
            length = self.length,
            sampleRate = self.sampleRate,
            stepRate = self.stepRate,
            randomState = pickle.dumps(self.rSource.get_state()),
            numpyVersion = N.__version__,
        )
        
        return AcqUI.cocoaCollection(d)
    
        
    
    
    def version(self):
        return 1
    
    @classmethod
    def stimulusID(self):
        return 'edu.washington.bwark.acqui.stim.stair_step_guassian'        
    


objc.removeAutoreleasePool()