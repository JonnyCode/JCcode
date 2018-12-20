classdef SplitNoiseStimulus < SpatialStimSuper
    properties %separate into private and public (params needed from IGOR)
        CenterShift     % number of pixels off center (+ shifts center right, - left)
        spotRadius      % radius of mask in pixels
        stimStd         % standard deviation of stim
        HighFreqCut     % highFrequency cuttoff
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
            
            % params
            CenterShift = params.CenterShift ;
            spotRadius = params.spotRadius ;  
            stimStd = params.stimStd ;
            HighFreqCut = params.HighFreqCut ;    
  
            % noise waveform 
            evalin('caller',['paramsStruct.randSeed=' num2str(params.randSeed) ';']); % setting random seed for recall later?
            randn('seed',params.randSeed) ;
            
            self.wave = normrnd(0,stimStd,1,self.spatial_stimPnts) ; % make vector for left side 
            self.wave = highPassFilter(self.wave,60,HighFreqCut) ; % high pass filter (signal,samplerate,frequency cutoff)
                
            % left and right side matricies
            self.lM = ones(self.screenY,(self.screenX/2)+CenterShift) *self.spatial_meanLevel ; 
            self.rM = ones(self.screenY,(self.screenX/2)-CenterShift) *self.spatial_meanLevel ;

            % mask
            if spotRadius~=0 ;              
                xmask = repmat([1:self.screenX],self.screenY,1)-(self.screenX/2) ;
                ymask = repmat([1:self.screenY]',1,self.screenX)-(self.screenY/2) ;
                self.maskPnts = sqrt(xmask.^2 + ymask.^2) > spotRadius;
            else 
                self.maskPnts = [] ;
            end
           
        end
               
        function drawParams = nextStimFrame(self)
            stimFrame = self.frameCounter - self.spatial_prepts + 1;
     
            M = [self.lM +self.wave(stimFrame),self.rM -self.wave(stimFrame)] ; % left side and right side have opposite and identical deviations around the mean 
            M(self.maskPnts) = self.spatial_meanLevel ; % mask as appropriate
            M_Scaled = gammaCorrect(M,self.gammaTable) ; % gamma correct matrix
            
            drawParams.tex = Screen('MakeTexture', self.windPtr, M_Scaled);

        end
    end
end



