classdef Layer < handle
    properties
        subunits
        input %stimulus or response of higher layer
        output %cell array representation suitabe for next layer or model response
        transform %class with "invert" and "revert" methods (to run on stimulus and filter)
        %layerType = []; %'input', 'middle', 'output' or combinations of the above
    end
    methods
        
        function run(self,options)
            L = length(self.subunits);
            self.output = cell(L,1);
            
            for i=1:L %for each subunit
                if ~isfield(options, 'preComputedInput')
                    self.subunits(i).input = self.input;
                end
                self.subunits(i).run(self.transform);
                self.output{i} = self.subunits(i).response;
            end
            
            if L==1 %if single response, take out of cell array
                self.output = self.output{1};
            end
        end
        
    end
end