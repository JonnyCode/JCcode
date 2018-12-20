classdef SpatialStimSuper < handle
    properties
        spatial_prepts
        spatial_stimpts 
        spatial_postpts 
        spatial_meanLevel        
        rand_behavior %0 = repeated seed, 1 = random seed, 2 = alternating
    end
    properties (Hidden = true)
        windPtr
        frameCounter = 1;
        texturePtr
        blankPtr
        screenX
        screenY
        callCounter
        randSeed
        initerror = 0;
    end
    
    methods
        function initSuper(self, windPtr, params)
            persistent callCounter;
            
            self.screenX = params.screenX;
            self.screenY = params.screenY;
            self.windPtr = windPtr;
            self.spatial_prepts = params.spatial_prepts;
            self.spatial_postpts = params.spatial_postpts;
            self.spatial_stimpts = params.spatial_stimpts;
            self.rand_behavior = params.rand_behavior;
            gammaTable = params.gammaTable;
                        
            if isfield(params,'spatial_meanLevel')
                self.spatial_meanLevel = gammaCorrect(params.spatial_meanLevel, gammaTable);
            else
                self.spatial_meanLevel = 0;
            end
            self.blankPtr = Screen('MakeTexture', self.windPtr, self.spatial_meanLevel);
            
            %deal with callCounter
            %disp(['InitStimulusSequence = ' num2str(params.InitStimulusSequence)]);
            if params.InitStimulusSequence %if init (first in loop) may change name of this param
                callCounter = 1;
            else
                callCounter = callCounter + 1;
            end
            self.callCounter = callCounter;
            disp(['CallCounter = ' num2str(self.callCounter)]);
            
            %set rand seed            
            if self.rand_behavior == 0 %repeat seed
                randSeed = 1;                 
            elseif self.rand_behavior == 1 %random seed
                randSeed = double(params.epochNumber);
            elseif self.rand_behavior == 2 %alternate
                if rem(self.callCounter,2)==1 %repeated
                    randSeed = 1;  
                else %random
                    randSeed = double(params.epochNumber+2); %add 2 so we don't get 0 or 1
                end
            end                
            self.randSeed = randSeed;
            evalin('caller',['params.randSeed=' num2str(randSeed) ';']);
            rand('seed', randSeed);
        end
        
        function drawParams = nextFrame(self)
            if self.frameCounter <= self.spatial_prepts; %in prepts
                drawParams.tex = self.blankPtr;
            elseif self.frameCounter > self.spatial_prepts+self.spatial_stimpts %in postpts
                drawParams.tex = self.blankPtr;
            elseif self.frameCounter > self.spatial_prepts; %in stimpts
                drawParams = self.nextStimFrame;
            end        
            self.frameCounter = self.frameCounter+1;
        end
        
    end
end

