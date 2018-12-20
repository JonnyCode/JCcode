classdef ProtocolClassSuper < handle
    properties
        AoNumber
        TrialDuration 
        NumRepeats
        TimeBetweenRepeats
        SampleRate
        TrialTag
    end
    properties (Hidden = true)
        OutputData
        LoopNumber
        AoNum
    end
    methods
        
        function self = SpatialStimulus(varargin)
         
        end
        

        
    end %methods
end %class