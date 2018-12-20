classdef ScanTextureStimulus < SpatialStimSuper
    properties %separate into private and public (params needed from IGOR)
        CorrFactor      % std of texture
        SpotRadius      % radius of mask
        ScanSpeed       % pixels/frame over which the image is scanned by the mask
        ScanAngle       % the angle over which the image is scanned by the mask
        RandTexture     % randomize the texture that is generated
        Contrast        % percent above and below .5 mean 
        HoldFrames      % frames between appearance of texture and movement
        SquareTexture   % 1 = make texure only black and white
    end
    properties (Hidden = true)
        Texture
        gammaTable
    end
    methods
        
        function self = ScanTextureStimulus(windPtr, params)   
            self.initSuper(windPtr, params) ;                
            
            pF = 1 ; % pixel factor: make texture smaller and expand (does not currently work properly!) 
            y = ceil(self.screenY/pF) ;
            x = ceil(self.screenX/pF) ;
            
            self.SpotRadius = params.SpotRadius/pF ;
            self.ScanSpeed = params.ScanSpeed/pF ;
            self.ScanAngle = params.ScanAngle ;
            self.Contrast = params.Contrast ;
            self.CorrFactor = params.CorrFactor/pF ;
            self.HoldFrames = params.HoldFrames ;
            
            ScanPts = self.spatial_stimpts*self.ScanSpeed + x - self.HoldFrames ; % number of frames per scan
            
            M = generateTexture2(ScanPts,y,self.CorrFactor,self.Contrast) ; % get texture
            
            if params.SquareTexture == 1 ;
                M(M>=.5)=0.5 + 0.5*self.Contrast ;
                M(M<.5)=0.5 - 0.5*self.Contrast ;
            end
            
            self.gammaTable = params.gammaTable ;
            M = gammaCorrect(M,self.gammaTable) ;
 
            % mask
            if self.SpotRadius~=0 ;              
                xmask = repmat([1:x],y,1)-(x/2) ;
                ymask = repmat([1:y]',1,x)-(y/2) ;
                maskPnts = sqrt(xmask.^2 + ymask.^2) > self.SpotRadius ;
            else 
                maskPnts = [] ;
            end    

            round=0 ;
            for a=1:self.ScanSpeed:ScanPts-x  ;
                round = round+1 ;
                Stim = M(:,a:a+x-1) ;
                Stim(maskPnts) = 0 ;
                self.Texture(round)=Screen('MakeTexture', self.windPtr, Stim) ;
            end
                
            [resident] = Screen('PreloadTextures', self.windPtr, self.Texture) ;
                     
        end
        
        function drawParams = nextStimFrame(self)
            frame = self.frameCounter-self.spatial_prepts ;
                                
            if frame<=self.HoldFrames ;
                drawParams.rotationAngle = self.ScanAngle ; % bar angle for this bar              
                drawParams.tex = self.Texture(1) ; 
            else
                drawParams.rotationAngle = self.ScanAngle ; % bar angle for this bar              
                drawParams.tex = self.Texture(frame-self.HoldFrames) ;
            end
            
        end
        
    end
end