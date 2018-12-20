classdef OdorPulseStimulus < ProtocolClassSuper
    properties %separate into private and public (user defined parameters)
        OdorPulseInitTime % (sec) start time of odor pulse
        OdorPulseDur % (sec) duration of odor pulse
        OdorPulseAmp % (volts) from daq to air valve
    end

    methods
        
        function OdorPulseStimulus
            self.OutputData = zeros(1,self.TrialDuration*self.SampleRate) ; 
            self.OutputData(self.OdorPulseInitTime*self.SampleRate:(self.OdorPulseInitTime+self.OdorPulseDur)*self.SampleRate,1)= self.OdorPulseAmp ;

               
        end  
        
    end
end