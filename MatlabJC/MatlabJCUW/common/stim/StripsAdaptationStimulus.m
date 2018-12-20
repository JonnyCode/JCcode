classdef StripsAdaptationStimulus < SpatialStimSuper
    properties
        stripWidth %pixels
        spotDiameter %pixels
        backgroundIntensity
        flashIntensity %can be a vector
        flashStartFrame %from beginning of stimpts
        flashDuration %frames
        flashPosition %0 = in phase w/bg, 1 = out of phase, 2 = alternating
        backgroundDiameter %pixels; 0 = full screen
    end
    properties (Hidden = true)
        flashIntensityList
        gammaTable
        Xpos %set to self.screenX/2        
        Ypos %set to self.screenY/2
        texturePtr_bg
        texturePtr_flash
    end
    
    methods
        
        function self = StripsAdaptationStimulus(windPtr, params)
            self.initSuper(windPtr, params);
            
            persistent flashIntensityList
                        
            self.stripWidth = params.stripWidth;
            self.spotDiameter = params.spotDiameter;
            self.flashStartFrame = params.flashStartFrame;
            self.flashDuration = params.flashDuration;            
            self.gammaTable = params.gammaTable;
            self.backgroundDiameter = params.backgroundDiameter;
            self.flashPosition = params.flashPosition;
            
            self.backgroundIntensity = gammaCorrect(params.backgroundIntensity,self.gammaTable);
            self.flashIntensity = gammaCorrect(params.flashIntensity,self.gammaTable);
            evalin('caller',['paramsStruct.flashIntensity_raw=[' num2str(self.flashIntensity) '];']);
            evalin('caller',['paramsStruct.backgroundIntensity_raw=[' num2str(self.backgroundIntensity) '];']);
                        
            %set Xpos and Ypos
            self.Xpos = round(self.screenX/2);
            self.Ypos = round(self.screenY/2);
            
            nItems = params.NumCycles;
            
            %multiple intensities             
            if length(self.flashIntensity) > 1
                if self.callCounter == 1 %if init
                    R = ceil(rand(1,nItems).*length(self.flashIntensity));
                    flashIntensityList = self.flashIntensity(R);
                    self.flashIntensityList = flashIntensityList;                    
                end
                %use the next value in the list
                self.flashIntensity = flashIntensityList(self.callCounter);
                evalin('caller',['paramsStruct.flashIntensity = ' num2str(self.flashIntensity) ';']); %save the actual intensity used
            end
            
            %make the texture
            spotRadius = round(self.spotDiameter/2);      
            bgRadius = round(self.backgroundDiameter/2); 
            
            %check radius (not checking bg radius)
            if self.Xpos-spotRadius <= 0 || self.Xpos+spotRadius > self.screenX || ...
                    self.Ypos-spotRadius <= 0 || self.Ypos+spotRadius > self.screenY
                disp('Spot Position error');
                self.initerror = 1;
                return;
            end
            
            %make the spot and save the textures
            
            %texture matrices 
            Mbg = ones(self.screenX,self.screenY).*self.spatial_meanLevel;
            Mflash = ones(self.screenX,self.screenY).*self.spatial_meanLevel;
            
            %set up the masks
            shortDim = min(self.screenX,self.screenY);
            square_bg = ones(shortDim,shortDim).*self.spatial_meanLevel;
            square_flash = square_bg;
            X = -shortDim/2+1:shortDim/2; 
            Y = X;
            Xmat = repmat(X,shortDim,1);
            Ymat = repmat(Y,shortDim,1)';
            if self.backgroundDiameter>0               
                circlePoints_bg = sqrt(Xmat.^2 + Ymat.^2)<bgRadius;  
            end
            circlePoints_flash = sqrt(Xmat.^2 + Ymat.^2)<spotRadius;            
            
            %add the bg strips to the bg square            
            bins = ceil(size(square_bg,1) ./ self.stripWidth);
            ind = 1;
            for b=1:bins
                %bg
                endInd = ind+self.stripWidth-1;
                if endInd>size(square_bg,1)
                    endInd = size(square_bg,1);
                end
                square_bg(ind:endInd,:) = self.backgroundIntensity;
                square_flash(ind:endInd,:) = self.backgroundIntensity;
                
                %flash                
                
                if self.flashPosition==0 || (self.flashPosition==2 && rem(self.callCounter,2) == 0)
                   %in phase w/bg 
                   startInd_flash = ind;
                   endInd_flash = endInd;
                else
                   %out of phase
                   startInd_flash = ind+self.stripWidth;
                   endInd_flash = startInd_flash+self.stripWidth-1;
                end
                
                if startInd_flash<size(square_bg,1)
                    if endInd_flash>size(square_bg,1)
                        endInd_flash = size(square_bg,1);
                    end
                    square_flash(startInd_flash:endInd_flash,:) = self.flashIntensity;
                end
                
                ind = ind+self.stripWidth*2;
            end            
            
            %mask off the circles
            if self.backgroundDiameter>0, square_bg(~circlePoints_bg) = self.spatial_meanLevel; end;
            square_flash(~circlePoints_flash) = square_bg(~circlePoints_flash);
            
            Mbg(self.Xpos-shortDim/2+1:self.Xpos+shortDim/2, self.Ypos-shortDim/2+1:self.Ypos+shortDim/2) = square_bg;
            Mflash(self.Xpos-shortDim/2+1:self.Xpos+shortDim/2, self.Ypos-shortDim/2+1:self.Ypos+shortDim/2) = square_flash;
            self.texturePtr_bg = Screen('MakeTexture', self.windPtr, Mbg');                        
            self.texturePtr_flash = Screen('MakeTexture', self.windPtr, Mflash');                        
        end
                
        function drawParams = nextStimFrame(self)
            if self.frameCounter - self.spatial_prepts >= self.flashStartFrame && ...
                    self.frameCounter - self.spatial_prepts < self.flashStartFrame + self.flashDuration
                drawParams.tex = self.texturePtr_flash;
            else
                drawParams.tex = self.texturePtr_bg;
            end
        end
        
    end
end
