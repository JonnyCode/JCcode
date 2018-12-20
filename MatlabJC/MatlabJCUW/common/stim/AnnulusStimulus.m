classdef AnnulusStimulus < SpatialStimSuper
    properties
        Xpos
        Ypos
        diameter %can be a list of diameters to choose from
        intensity
        integratedIntensity %total intensity of spot, 0 to ignore and use intensity instead
    end
    properties (Hidden = true)
        %persistent properties
        diameterList
        XposList
        YposList
        intensityList
        integratedIntensityList
        gammaTable
    end
    
    methods
        
        function self = AnnulusStimulus(windPtr, params)
            self.initSuper(windPtr, params);
            
            self.Xpos = params.Xpos ;
            self.Ypos = params.Ypos ;
            self.diameter = params.diameter ;
            self.intensity = gammaCorrect(params.intensity,self.gammaTable) ;
            evalin('caller',['paramsStruct.intensity_raw=[' num2str(self.intensity) '];']);
            
            self.intergratedIntensity = params.integratedIntensity ;
                     
            %multiple diameters
            if length(self.diameter) > 1
                if self.callCounter == 1 %if init
                    self.diameterList = shuffle(self.diameter);
                end
                %use the next value in the list
                self.diameter = self.diameterList(self.callCounter);
                %params.diameter = self.diameter;
                evalin('caller',['paramsStruct.diameter = ' num2str(self.diameter) ';']); %save the actual diameter used
            end
            
            %multiple X,Y positions
            if length(self.Xpos) > 1
                if self.callCounter == 1 %if init
                    self.XposList = shuffle(self.Xpos);
                    self.YposList = shuffle(self.Ypos); 
                end
                %use the next value in the list
                self.Xpos = XposList(self.callCounter);
                self.Ypos = YposList(self.callCounter);

                evalin('caller',['paramsStruct.Xpos = ' num2str(self.Xpos) ';']); %save the actual diameter used
                evalin('caller',['paramsStruct.Ypos = ' num2str(self.Ypos) ';']); %save the actual diameter used
            end
            
            %multiple intensities 
            if length(self.intensity) > 1
                if self.callCounter == 1 %if init
                    self.intensityList = shuffle(self.intensity);
                end
                %use the next value in the list
                self.intensity = self.intensityList(self.callCounter);
                %params.diameter = self.diameter;
                evalin('caller',['paramsStruct.intensity = ' num2str(self.intensity) ';']); %save the actual diameter used
            end
            
            %multiple integrated intensities
            if length(self.integratedIntensity) > 1
                if self.callCounter == 1 %if init
                    self.integratedIntensityList = shuffle(self.integratedIntensity);
                end
                %use the next value in the list
                self.integratedIntensity = self.integratedIntensityList(self.callCounter);
                %params.diameter = self.diameter;
                evalin('caller',['paramsStruct.integratedIntensity = ' num2str(self.integratedIntensity) ';']); %save the actual diameter used
            end
            
            
            %make the texture
            radius = round(self.diameter/2);
            X = -radius:radius;
            Y = -radius:radius;
            Xmat = repmat(X,2*radius+1,1);
            Ymat = repmat(Y,2*radius+1,1)';
            circlePoints = sqrt(Xmat.^2 + Ymat.^2)<radius;
            
            if ~self.integratedIntensity
                M = ones(self.screenX,self.screenY).*self.intensity;
                square = ones(2*radius+1,2*radius+1).*self.intensity;
                square(circlePoints) = self.spatial_meanLevel;
            else
                nPixels = sum(sum(circlePoints));
                Ivalue = self.integratedIntensity/nPixels;
                disp(['Diamteter: ' num2str(self.diameter) ' Ivalue: ' num2str(Ivalue)]);
                if Ivalue < 0 || Ivalue > 1
                    error(['Ivalue of ' num2str(Ivalue) 'is out of range']);
                    self.initerror = 1;
                    return;
                end
                pixelValue = gammaCorrect(Ivalue,self.gammaTable);
                disp(['PixelValue: ' num2str(pixelValue)]);
                evalin('caller',['paramsStruct.intensity_raw=' num2str(pixelValue) ';']);
                M = ones(self.screenX,self.screenY).*pixelValue;
                square = ones(2*radius+1,2*radius+1).*pixelValue;
                square(circlePoints) = self.spatial_meanLevel;
            end
            
            M(self.Xpos-radius:self.Xpos+radius, self.Ypos-radius:self.Ypos+radius) = square;
            self.texturePtr = Screen('MakeTexture', self.windPtr, M');          
            
        end
        
        function drawParams = nextStimFrame(self)
            drawParams.tex = self.texturePtr;
        end
        
    end
end