classdef CheckerboardGreyStimulus < SpatialStimSuper
    properties %separate into private and public (params needed from IGOR)
        checkSizeX
        checkSizeY
        frameDwell %frames for each check pattern
        IntensityMean
        IntensityStd
        nCheckMaskX 
        nCheckMaskY
    end
    properties (Hidden = true)
        nChecksX
        nChecksY
        randStream
        gammaTable
        currentTexture
    end
    methods
        
        function self = CheckerboardGreyStimulus(windPtr, params)

            self.initSuper(windPtr, params);
            
            self.checkSizeX = params.checkSizeX;
            self.checkSizeY = params.checkSizeY;
            self.frameDwell = params.frameDwell;
            self.IntensityMean = params.IntensityMean ;
            self.IntensityStd = params.IntensityStd ;
            self.nCheckMaskX = params.nCheckMaskX ;
            self.nCheckMaskY = params.nCheckMaskY ;
            
            self.gammaTable = params.gammaTable ;
            
            self.nChecksX = self.screenX./ self.checkSizeX ;
            self.nChecksY = self.screenY./ self.checkSizeY ;
            
            evalin('caller',['paramsStruct.randSeed=' num2str(params.randSeed) ';']);
            randn('seed',params.randSeed) ;
            
        end
        
        function drawParams = nextStimFrame(self)
            stimFrame = self.frameCounter - self.spatial_prepts - 1;
            if rem(stimFrame, self.frameDwell) == 0
                
                %make new pattern
                M = normrnd(self.IntensityMean,self.IntensityStd,self.nChecksY,self.nChecksX) ;
                M_masked = M ;
                M_masked([1:self.nChecksMaskY,end-self.nChecksMaskY+1:end],[1:self.nChecksMaskX,end-self.nChecksMaskX+1:end]) = self.IntensityMean ;
                
                M_Scaled = gammaCorrect(M_masked,self.gammaTable) ;  
                
                self.currentTexture = Screen('MakeTexture', self.windPtr, M_Scaled);
            end
            drawParams.tex = self.currentTexture;
            %drawParams.rotationAngle = self.frameCounter;
        end
        
    end
end