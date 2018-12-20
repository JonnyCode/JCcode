classdef SplitFieldStimulus < SpatialStimSuper
    properties %separate into private and public (params needed from IGOR)
        intensityRight    % Bar cannot be wider than 424 pixels
        intensityLeft   % Bar cannot be larger than 212 pixels
        preFramesRight    % pixels/frame
        preFramesLeft % number of frames between next bar angle presented
<<<<<<< .mine
        spotRadius
=======
        spotRadius % number of pixels before mask (diameter of split field)
>>>>>>> .r544
    end
    properties (Hidden = true)
        i
        M
        MR
        ML
        MRL
        textureM
        textureMR
        textureML
        textureMRL

    end
    methods
        
        function self = SplitFieldStimulus(windPtr, params)
            self.initSuper(windPtr, params) ;
            
            self.intensityRight = gammaCorrect(params.intensityRight,params.gammaTable);
            self.intensityLeft = gammaCorrect(params.intensityLeft,params.gammaTable);
            
            self.preFramesRight = params.preFramesRight ;
            self.preFramesLeft = params.preFramesLeft ;
            
            self.spotRadius = params.spotRadius ;
            
            self.i = rem(self.callCounter,length(self.intensityRight)) ;
            if self.i==0 ;
                self.i=length(self.intensityRight);
            end
            
            if self.rand_behavior==1 ;
                R = randperm(length(self.intensityRight));
                self.i = R(self.i) ;
            end
            
            evalin('caller',['paramsStruct.intensityRight = ' num2str(self.intensityRight(self.i)) ';']); %save the actual diameter used
            evalin('caller',['paramsStruct.intensityLeft = ' num2str(self.intensityLeft(self.i)) ';']); %save the actual diameter used 
            
            % stim matricies
            self.M = ones(self.screenY,self.screenX).*self.spatial_meanLevel;
            self.MR = self.M ;
            self.ML = self.M ;
            self.MRL = self.M ;
            
            self.MR(:,self.screenX/2+1:end) = self.intensityRight(self.i) ;
            self.ML(:,1:self.screenX/2) = self.intensityLeft(self.i) ;
            self.MRL(:,self.screenX/2+1:end) = self.intensityRight(self.i) ;
            self.MRL(:,1:self.screenX/2) = self.intensityLeft(self.i) ;
            
            % mask
            if self.spotRadius~=0 ;
                
                xmask = repmat([1:self.screenX],self.screenY,1)-(self.screenX/2) ;
                ymask = repmat([1:self.screenY]',1,self.screenX)-(self.screenY/2) ;
                maskPnts = sqrt(xmask.^2 + ymask.^2) > self.spotRadius(self.i) ;

                self.MR(maskPnts) = self.spatial_meanLevel ;
                self.ML(maskPnts) = self.spatial_meanLevel ; 
                self.MRL(maskPnts) = self.spatial_meanLevel ;
            end
            
            self.textureM = Screen('MakeTexture', self.windPtr, self.M) ;
            self.textureMR = Screen('MakeTexture', self.windPtr, self.MR) ;
            self.textureML = Screen('MakeTexture', self.windPtr, self.ML) ;
            self.textureMRL = Screen('MakeTexture', self.windPtr, self.MRL) ;
            
        end
        
        function drawParams = nextStimFrame(self)
            frame = self.frameCounter-self.spatial_prepts ;
            texture = self.textureM ;
            
            if frame>self.preFramesRight(self.i) ;
                texture = self.textureMR ;
            end
            
            if frame>self.preFramesLeft(self.i) ;
                texture = self.textureML ;
            end
            
            if frame>self.preFramesRight(self.i) & frame>self.preFramesLeft(self.i) ;
                texture = self.textureMRL ;
            end
            
            drawParams.tex = texture;
        end
    end
end
