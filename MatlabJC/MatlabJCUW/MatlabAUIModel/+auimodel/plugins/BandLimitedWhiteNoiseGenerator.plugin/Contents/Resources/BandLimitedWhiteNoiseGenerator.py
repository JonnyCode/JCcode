#
#  BandLimitedWhiteNoiseGenerator.py
#  AcqUI
#
#  Created by Barry Wark on 2/22/06.
#  Copyright (c) 2006 Barry Wark. All rights reserved.
#

from Foundation import NSLog, NSDictionary
from acqui.AcqUI import AUIStimulusPyBase

from acqui.noise_source import BandLimitedNoiseSource, GaussianNoiseSource
from acqui import AcqUI

import cPickle as pickle
import numpy

import objc

BLWN_MEAN_KEY = 'mean'
BLWN_VAR_KEY = 'var'
BLWN_FREQ_KEY = 'fCutoff'
BLWN_NPOLES_KEY = 'nPoles'
BLWN_SAMPLE_RATE_KEY = 'sampleRate'
BLWN_LENGTH_KEY = 'length'
BLWN_RANDOM_STATE_KEY = 'rState'

class BandLimitedWhiteNoiseGenerator (AUIStimulusPyBase):
    
    
    def init(self):
        self = super(BandLimitedWhiteNoiseGenerator, self).init()
        
        if(self != None):
            self.mean = 0
            self.var = 1
            self.nPoles = 1
            self.fCutoff = 1000
            self.length = 1
            self.sampleRate = 1e4
        
        return self
    
    def stimulusDataForParamsDictionary_version_(self, p, v):
        mean = p[BLWN_MEAN_KEY].doubleValue()
        var = p[BLWN_VAR_KEY].doubleValue()
        len = p[BLWN_LENGTH_KEY].doubleValue()
        hz = p[BLWN_SAMPLE_RATE_KEY].doubleValue()
        fCutoff = p[BLWN_FREQ_KEY].doubleValue()
        nPoles = p[BLWN_NPOLES_KEY].intValue()
        
        if(numpy.__version__ != p['numpyVersion']):
            NSLog('Warning: A different numpy version (' + p['numpyVersion'] + ') was used to create BLWN stimulus. Current version: ' + numpy.__version__)
        
        self.rSource.set_state(pickle.loads(str(p[BLWN_RANDOM_STATE_KEY])))
        
        ns = BandLimitedNoiseSource(GaussianNoiseSource(mean, var, self.rSource), fCutoff, nPoles)
        
        if(v==1): #random phase construction
            arr = ns.noise(len, hz, rPhase=True)
        else:
            arr = ns.noise(len, hz)
        
        return arr.astype(AcqUI.dtypeForAnalogOutput())
    
    
    def paramsDictionary(self):
        self.fCutoff = min(self.fCutoff, self.sampleRate/2.)
        
        return NSDictionary.dictionaryWithObjectsAndKeys_(
                    self.mean, BLWN_MEAN_KEY, 
                    self.var, BLWN_VAR_KEY, 
                    self.fCutoff, BLWN_FREQ_KEY, 
                    self.nPoles, BLWN_NPOLES_KEY, 
                    self.sampleRate, BLWN_SAMPLE_RATE_KEY, 
                    self.length, BLWN_LENGTH_KEY, 
                    pickle.dumps(self.rSource.get_state()[:3]), BLWN_RANDOM_STATE_KEY, 
                    numpy.__version__, 'numpyVersion', 
                    None)
    
    
    def version(self):
        return 2
        
    @classmethod
    def stimulusID(self):
        return 'edu.washington.bwark.acqui.stim.BandLimitedWhiteNoise'
    


objc.removeAutoreleasePool()
