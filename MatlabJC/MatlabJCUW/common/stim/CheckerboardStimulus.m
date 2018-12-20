classdef CheckerboardStimulus < SpatialStimSuper
    properties %separate into private and public (params needed from IGOR)
        checkSizeX
        checkSizeY
        frameDwell %frames for each check pattern
        intensity
    end
    properties (Hidden = true)
        nChecksX
        nChecksY
        randStream
        gammaTable
        currentTexture
    end
    methods
        
        function self = CheckerboardStimulus(windPtr, params)
            self.initSuper(windPtr, params);
            
            self.checkSizeX = params.checkSizeX;
            self.checkSizeY = params.checkSizeY;
            self.frameDwell = params.frameDwell;
            self.gammaTable = params.gammaTable;
            self.intensity = gammaCorrect(params.intensity,self.gammaTable);
            evalin('caller',['paramsStruct.intensity_raw=' num2str(self.intensity) ';']);
            evalin('caller',['paramsStruct.randSeed=' num2str(params.randSeed) ';']);
            
            self.nChecksX = self.screenX./ self.checkSizeX;
            self.nChecksY = self.screenY./ self.checkSizeY;
        end
        
        function drawParams = nextStimFrame(self)
            stimFrame = self.frameCounter - self.spatial_prepts - 1;
            if rem(stimFrame, self.frameDwell) == 0
                %make new pattern
                M = rand(self.nChecksY,self.nChecksX)>.5;
                M_scaled = double(M)*self.intensity;
                self.currentTexture = Screen('MakeTexture', self.windPtr, M_scaled);
            end
            drawParams.tex = self.currentTexture;
            %drawParams.rotationAngle = self.frameCounter;
        end
        
    end
end