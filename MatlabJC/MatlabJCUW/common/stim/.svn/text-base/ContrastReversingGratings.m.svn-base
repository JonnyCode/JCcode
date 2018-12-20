classdef ContrastReversingGratings < SpatialStimSuper
    properties
        gratingWidth %pixels in one period
        spotDiameter %pixels
        amplitude %can be a vector: stim range is meanLevel +- amplitide
        stimPeriod %frames, must be even
    end
    properties (Hidden = true)
        amplitudeList
        gratingWidthList
        gammaTable
        Xpos %set to self.screenX/2        
        Ypos %set to self.screenY/2
        texturePtrVec
    end
    
    methods
        
        function self = ContrastReversingGratings(windPtr, params)
            self.initSuper(windPtr, params);
            
            persistent amplitudeList
            persistent gratingWidthList
                        
            self.gratingWidth = params.gratingWidth;
            self.spotDiameter = params.spotDiameter;
            self.gammaTable = params.gammaTable;
            self.stimPeriod = params.stimPeriod;
            
            self.amplitude = gammaCorrect(params.amplitude,self.gammaTable);
            evalin('caller',['paramsStruct.amplitude_raw=[' num2str(self.amplitude) '];']);
                        
            %set Xpos and Ypos
            self.Xpos = round(self.screenX/2);
            self.Ypos = round(self.screenY/2);
            
            nItems = params.NumCycles;
            
            %multiple amplitudes             
            if length(self.amplitude) > 1
                if self.callCounter == 1 %if init
                    R = ceil(rand(1,nItems).*length(self.amplitude));
                    amplitudeList = self.amplitude(R);
                    self.amplitudeList = amplitudeList;                    
                end
                %use the next value in the list
                self.amplitude = amplitudeList(self.callCounter);
                evalin('caller',['paramsStruct.amplitude = ' num2str(self.amplitude) ';']); %save the actual amplitude used
            end
            
            %multiple grating widths
            if length(self.gratingWidth) > 1
                if self.callCounter == 1 %if init
                    R = ceil(rand(1,nItems).*length(self.gratingWidth));
                    gratingWidthList = self.gratingWidth(R);
                    self.gratingWidthList = gratingWidthList;                    
                end
                %use the next value in the list
                self.gratingWidth = gratingWidthList(self.callCounter);
                evalin('caller',['paramsStruct.gratingWidth = ' num2str(self.gratingWidth) ';']); %save the actual grating width used
            end
                        
            %make the texture
            spotRadius = round(self.spotDiameter/2);      
            
            %check radius (not checking bg radius)
            if self.Xpos-spotRadius <= 0 || self.Xpos+spotRadius > self.screenX || ...
                    self.Ypos-spotRadius <= 0 || self.Ypos+spotRadius > self.screenY
                disp('Spot Position error');
                self.initerror = 1;
                return;
            end
            
            
            %for each frame in the period
            for f=1:self.stimPeriod/2
                %make the spot and save the textures
                
                %texture matrices
                Mflash = ones(self.screenX,self.screenY).*self.spatial_meanLevel;
                
                %set up the masks
                shortDim = min(self.screenX,self.screenY);
                square_flash = ones(shortDim,shortDim).*self.spatial_meanLevel;
                X = -shortDim/2+1:shortDim/2;
                Y = X;
                Xmat = repmat(X,shortDim,1);
                Ymat = repmat(Y,shortDim,1)';
                circlePoints_flash = sqrt(Xmat.^2 + Ymat.^2)<spotRadius;
                
                %add the sinewave grating
                c1Amp = (f-1)/(self.stimPeriod/2-1);
                c2Amp = 1-c1Amp;
                
                cycles = ceil(size(square_flash,1) ./ self.gratingWidth);
                square_flash_temp = self.spatial_meanLevel + repmat(self.amplitude.*...
                    (c1Amp*sin(linspace(0,2*pi,self.gratingWidth)) - c2Amp*sin(linspace(0,2*pi,self.gratingWidth))),...
                    size(square_flash,2),cycles)';
                square_flash = square_flash_temp(1:size(square_flash,1),1:size(square_flash,2));
                
                %mask off the circles
                square_flash(~circlePoints_flash) = self.spatial_meanLevel;
                
                Mflash(self.Xpos-shortDim/2+1:self.Xpos+shortDim/2, self.Ypos-shortDim/2+1:self.Ypos+shortDim/2) = square_flash;
                self.texturePtrVec(f) = Screen('MakeTexture', self.windPtr, Mflash');
                
                
            end
        end
                
        function drawParams = nextStimFrame(self)
            R = rem(self.frameCounter - self.spatial_prepts, self.stimPeriod);
            if R<self.stimPeriod/2
                phaseInd = R+1;
            else %reverse
                phaseInd = self.stimPeriod - R;
            end
            drawParams.tex = self.texturePtrVec(phaseInd);
        end
    end
end
