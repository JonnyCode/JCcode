#!/usr/bin/env python
# encoding: utf-8
"""
ConcatenateStimulusGenerator.py

Created by  on 2006-08-03.
Copyright (c) 2006 Barry Wark. All rights reserved.
"""

from acqui.AcqUI import AUICompoundStimulusGenerator

import numpy as N

import objc


class ConcatenateStimulusGenerator (AUICompoundStimulusGenerator):
    
    def getOrderedStimuli(self):
        """getOrderedStim"""
        assert(N.iterable(self.subStimuli))
        return self.subStimuli
    
    def setOrderedStimuli(self, stim):
        """setOrderedStim"""
        assert(N.iterable(stim))
        self.subStimuli = list(stim)
    
    orderedStimuli = property(fget=getOrderedStimuli, fset=setOrderedStimuli)
    
    def init(self):
        self = super(ConcatenateStimulusGenerator, self).init()
        
        if(self != None):
            self.combFn = lambda (s): N.hstack(s)
            self.lengthFn = sum
        
        return self
    
    
    def version(self):
        return 2
    
    
    @classmethod
    def stimulusID(self):
        return 'edu.washington.bwark.acqui.stim.concatenate'        
    

objc.removeAutoreleasePool()
