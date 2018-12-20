classdef MovingBarStimulus < SpatialStimSuper
    properties %separate into private and public (params needed from IGOR)
        BarWidth    % Bar cannot be wider than 424 pixels
        BarHeight   % Bar cannot be larger than 212 pixels
        BarSpeed    % pixels/frame
        InterBarPts % number of frames between next bar angle presented
        BarAngle    % vector of desired bar angles
        randBarAngle % 1 = randomize the order of input bar angles
        intensity
        executableString %for Igor to sepcify stimpts
    end
    properties (Hidden = true)
        randStream
        gammaTable
        BarPts
        BarTexture
        BlankTexture
        maskLength
    end
    methods
        
        function self = MovingBarStimulus(windPtr, params)
            %could do this automatically (set each prop from params in
            %loop)
            pF = 8 ;
            
            self.initSuper(windPtr, params) ;                
            
            self.BarWidth = floor(params.BarWidth/pF) ;
            self.BarHeight = floor(params.BarHeight/pF) ;
            self.BarSpeed = floor(params.BarSpeed/pF) ;
            self.InterBarPts = params.InterBarPts ;
            
            self.maskLength = floor(2*(sqrt(((min(params.screenY/pF,params.screenX/pF)/2)^2) - (((self.BarHeight)/2)^2)))) ; % length of rectangular mask through which the bar will presented (ie distance bar travels)  
            
            if params.randBarAngle == 1 ; % if you want to use the input BarAngles in a random order 
                params.BarAngle = Shuffle(params.BarAngle)' ;  
                evalin('caller', ['paramsStruct.BarAngle = [' num2str(params.BarAngle') '];']) ; % sends out to paramsStruct
            end
            self.BarAngle = params.BarAngle ;                
            
            self.BarPts = ceil((self.maskLength+self.BarWidth-1)/self.BarSpeed) ; % number of frames per bar 
                  
            self.gammaTable = params.gammaTable ;
            self.intensity = gammaCorrect(params.intensity,self.gammaTable) ;
            evalin('caller',['paramsStruct.intensity_raw=' num2str(self.intensity) ';']) ;
            
            TotalStmPts = (self.BarPts+self.InterBarPts)*length(self.BarAngle) ; % number of stim points is determined by bar speed, number of angles, and inter bar pause
            if self.spatial_stimpts ~= TotalStmPts ;
                msgbox(['matlab stm pnts =',num2str(TotalStmPts),'and Igor stm pnts =',num2str(self.spatial_stimpts)])
            end
            self.spatial_stimpts = TotalStmPts;
            evalin('caller', ['paramsStruct.spatial_stimpts =' num2str(TotalStmPts) ';']) ; % sends out to paramsStruct

            centerYFactor = floor((params.screenY/pF-self.BarHeight)/2) ;  
            centerXFactor = floor((params.screenX/pF-self.maskLength)/2) ; 
            
            startPnts = [ones(1,self.BarWidth),[2:self.maskLength]] + centerXFactor ; 
            endPnts = [[1:self.maskLength],ones(1,self.BarWidth-1)*self.maskLength] + centerXFactor ;
           
            Blank = zeros(params.screenY/pF,params.screenX/pF) ;
          
            for a=1:self.BarPts ;
                Bar = Blank ;
                Bar(centerYFactor:centerYFactor+self.BarHeight, startPnts(a*self.BarSpeed-self.BarSpeed+1):endPnts(a*self.BarSpeed-self.BarSpeed+1)) = self.intensity ; % make matrix of bar moving
                self.BarTexture(a)=Screen('MakeTexture', self.windPtr, Bar) ;
            end
                self.BlankTexture = Screen('MakeTexture',self.windPtr, Blank) ;
                
                [resident] = Screen('PreloadTextures', self.windPtr, self.BarTexture)
                     
        end
        
        function drawParams = nextStimFrame(self)
            frame = self.frameCounter-self.spatial_prepts ;
                        
            if rem(frame,(self.BarPts+self.InterBarPts))>self.BarPts || rem(frame,(self.BarPts+self.InterBarPts))==0 ; %if your during a inter bar point...
                
                drawParams.tex = self.BlankTexture ; 
            else
                           
                drawParams.rotationAngle = self.BarAngle(ceil(frame/(self.BarPts+self.InterBarPts))) ; % bar angle for this bar
                
                drawParams.tex = self.BarTexture(rem(frame,(self.BarPts+self.InterBarPts))) ;
                
            end
            %imagesc(Bar);
            %pause;
                  
            
        end
        
    end
end