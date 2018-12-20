classdef SplitNoiseStimulus < SpatialStimSuper
    properties %separate into private and public (params needed from IGOR)
        CenterShift     % number of pixels off center (+ shifts center right, - left)
        spotRadius      % radius of mask
        stimStd         % standard deviation of stim
        FrameDwell      % number of frames 
    end
    properties (Hidden = true)
        wave
        rM
        lM
        maskPnts

    end
    methods
        
        function self = SpatialGradingStimulus(windPtr, params)
            self.initSuper(windPtr, params) ;
            
            CenterShift = params.CenterShift ;
            spotRadius = params.spotRadius ;  
            stimStd = params.stimStd ;
            self.FrameDwell = params.FrameDwell ;    
  
            % waveform on left side
            evalin('caller',['paramsStruct.randSeed=' num2str(params.randSeed) ';']); % setting random seed for recall later?
            randn('seed',params.randSeed) ;
            
            self.wave = normrnd(self.spatial_meanLevel,stimStd,1,self.spatial_stimPnts) ; % make vector for left side
            self.wave = gammaCorrect(self.wave,params.gammaTable) ; % gamma correct wave    
                
            % left and right side matricies
            lM = ones(self.screenY,(self.screenX/2)+CenterShift) ; 
            rM = ones(self.screenY,(self.screenX/2)-CenterShift) ;
            
            % mask
            if spotRadius(self.index)~=0 ;              
                xmask = repmat([1:self.screenX],self.screenY,1)-(self.screenX/2) ;
                ymask = repmat([1:self.screenY]',1,self.screenX)-(self.screenY/2) ;
                self.maskPnts = sqrt(xmask.^2 + ymask.^2) > spotRadius(self.index) ;
            else 
                self.maskPnts = [] ;
            end
           
        end
       
        
        function drawParams = nextStimFrame(self)
            stimFrame = self.frameCounter - self.spatial_prepts - 1;
            
            if rem(stimFrame, self.frameDwell) == 0 % make a new pattern if your not dwelling
                
                M = [lM*wave(stimFrame),-rM*wave(stimFrame)] ; % left side and right side are identical but right side is negative  
                M(maskPnts) = self.spatial_meanLevel ; % mask as appropriate
                
                self.currentTexture = Screen('MakeTexture', self.windPtr, M);
            
            end
            drawParams.tex = self.currentTexture;

        end
    end
end



