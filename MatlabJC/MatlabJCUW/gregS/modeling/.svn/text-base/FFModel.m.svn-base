classdef FFModel < handle %feedforward model
    properties
        layers %first layer is input layer, last is output layer
        %connectivityRule %function handle, takes Model as parameter
        output %cell array of output for each "trial". Trials are identical except for noise
        nTrials
        stimulus %in whatever representation it comes (already transformed)
        options %structure with optional parameters (that affect run method)
    end
    methods
        function run(self)
            %self.connectivityRule(self); %set up network connectivity
            
            L = length(self.layers);
            self.output = cell(self.nTrials,1);
            
            for i=1:self.nTrials
                disp(['Running Model: Trial ' num2str(i) ' of ' num2str(self.nTrials)]);
                self.layers(1).input = self.stimulus;
                
                %run each layer
                for j=1:L
                    self.layers(j).run(self.options);
                    if j<L
                        self.layers(j+1).input = self.layers(j).output;
                    end
                end
                
                %collect output
                self.output{i} = self.layers(end).output;
            end
            if self.nTrials == 1 %get rid of cell array container
                self.output = self.output{i};
            end
        end
    end
end