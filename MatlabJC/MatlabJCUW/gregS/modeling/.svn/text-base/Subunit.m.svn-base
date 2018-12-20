classdef Subunit < handle
    properties
        %if top level unit, input is stimulus, so we need...
        filter %normalized gain for each spatial and temporal "dimension". This is a vector, already transformed.
        
        %optional: only for level 1 units for spatial pre-filtering
        spatialPreFilter

        %if not top level unit
        inputGain %vector of normalized gains (can be zero) from higher level
                
        %any unit
        input %stim or responses from previous level
        outputNL %function handle (some default for linear), can be a function of which output (maybe?) or not
        response %gets set when model is run (hopefully scalar or vector)
    end
    methods
        
        function run(self, transform)
            %transform step
            
            if isempty(self.inputGain) %top (stimulus) level unit                
                %stimulus comes pre-transformed                                
                self.response = self.outputNL(self.filter.*self.input, transform);                
            else
                %sum over inputs
                for i=1:length(self.input)
                    if i==1
                        R = self.input{i}.*self.inputGain(i);
                    else
                        R = R + self.input{i}.*self.inputGain(i);
                    end
                end
                self.response = self.outputNL(self.filter.*R, transform);
            end
        end       
    end
end