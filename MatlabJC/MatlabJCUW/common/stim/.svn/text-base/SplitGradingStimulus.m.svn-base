classdef SplitGradingStimulus < SpatialStimSuper
    properties %separate into private and public (params needed from IGOR)
        BarWidth    % pixel width of half cycle of square wave
        Phase    % fraction of cycle to move square wave over
        intensityRight    % 
        intensityLeft   % 
        preFramesRight    % 
        preFramesLeft % 
        spotRadius
    end
    properties (Hidden = true)
        i
        textureM
        textureMR
        textureML
        textureMRL

    end
    methods
        
        function self = SplitGradingStimulus(windPtr, params)
            self.initSuper(windPtr, params) ;
            
            self.intensityRight = gammaCorrect(params.intensityRight,params.gammaTable);
            self.intensityLeft = gammaCorrect(params.intensityLeft,params.gammaTable);
            
            self.preFramesRight = params.preFramesRight ;
            self.preFramesLeft = params.preFramesLeft ;
            
            self.BarWidth = params.BarWidth ;
            self.Phase = params.Phase ;
            self.spotRadius = params.spotRadius ;
           
            self.i = rem(self.callCounter,length(self.intensityRight)) ;
            if self.i==0 ;
                self.i=length(self.intensityRight);
            end
            
            sineWave = sin((pi/self.BarWidth(self.i))*[1:self.screenX]+self.Phase(self.i)*pi/180) ;
            sineWave = repmat(sineWave,self.screenY,1) ;
                     
            if self.rand_behavior==1 ;
                R = randperm(length(self.intensityRight));
                self.i = R(self.i) ;
            end
            
            M = ones(self.screenY,self.screenX).*self.spatial_meanLevel;
            MR = M ;
            ML = M ;
            MRL = M ;
            
            MR(sineWave>0)= self.intensityRight(self.i) ;
            ML(sineWave<0)= self.intensityLeft(self.i) ;
            
            MRL(sineWave>0)= self.intensityRight(self.i) ;
            MRL(sineWave<0)= self.intensityLeft(self.i) ;
            
            % mask
            if self.spotRadius~=0 ;
                
                xmask = repmat([1:self.screenX],self.screenY,1)-(self.screenX/2) ;
                ymask = repmat([1:self.screenY]',1,self.screenX)-(self.screenY/2) ;
                maskPnts = sqrt(xmask.^2 + ymask.^2) > self.spotRadius(self.i) ;

                MR(maskPnts) = self.spatial_meanLevel ;
                ML(maskPnts) = self.spatial_meanLevel ; 
                MRL(maskPnts) = self.spatial_meanLevel ;
            end
            
            self.textureM = Screen('MakeTexture', self.windPtr, M) ;
            self.textureMR = Screen('MakeTexture', self.windPtr, MR) ;
            self.textureML = Screen('MakeTexture', self.windPtr, ML) ;
            self.textureMRL = Screen('MakeTexture', self.windPtr, MRL) ;
            
            
            evalin('caller',['paramsStruct.intensityRight = ' num2str(self.intensityRight(self.i)) ';']); %save the actual diameter used
            evalin('caller',['paramsStruct.intensityLeft = ' num2str(self.intensityLeft(self.i)) ';']); %save the actual diameter used            
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