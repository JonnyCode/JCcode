classdef SingleSpotStimulus < SpatialStimSuper
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
        
        function self = SingleSpotStimulus(windPtr, params)
            self.initSuper(windPtr, params);
            
            persistent diameterList
            persistent intensityList
            persistent integratedIntensityList
            persistent XposList
            persistent YposList
            
            self.Xpos = params.Xpos;
            self.Ypos = params.Ypos;
            self.diameter = params.diameter;
            self.gammaTable = params.gammaTable;
            self.intensity = gammaCorrect(params.intensity,self.gammaTable);
            evalin('caller',['paramsStruct.intensity_raw=[' num2str(self.intensity) '];']);
            self.integratedIntensity = params.integratedIntensity;
            
            nItems = params.NumCycles;
            
            %multiple diameters
            if length(self.diameter) > 1
                if self.callCounter == 1 %if init
                    R = ceil(rand(1,nItems).*length(self.diameter));
                    diameterList = self.diameter(R);
                    self.diameterList = diameterList;
                end
                %use the next value in the list
                self.diameter = diameterList(self.callCounter);
                %params.diameter = self.diameter;
                evalin('caller',['paramsStruct.diameter = ' num2str(self.diameter) ';']); %save the actual diameter used
            end
            
            %multiple X,Y positions
            if length(self.Xpos) > 1
                if self.callCounter == 1 %if init
                    if length(self.Xpos) ~= length(self.Ypos)
                        error('Xpos and Ypos must be the same length');
                        self.initerror = 1;
                        return;
                    end
                    R = ceil(rand(1,nItems).*length(self.Xpos));
                    XposList = self.Xpos(R);
                    YposList = self.Ypos(R);
                    self.XposList = XposList;
                    self.YposList = YposList;
                end
                %use the next value in the list
                self.Xpos = XposList(self.callCounter);
                params.Xpos = self.Xpos;
                self.Ypos = YposList(self.callCounter);
                params.Ypos = self.Ypos;
                evalin('caller',['paramsStruct.Xpos = ' num2str(self.Xpos) ';']); %save the actual diameter used
                evalin('caller',['paramsStruct.Ypos = ' num2str(self.Ypos) ';']); %save the actual diameter used
            end
            
            %multiple intensities 
            if length(self.intensity) > 1
                if self.callCounter == 1 %if init
                    R = ceil(rand(1,nItems).*length(self.intensity));
                    intensityList = self.intensity(R);
                    self.intensityList = intensityList;
                end
                %use the next value in the list
                self.intensity = intensityList(self.callCounter);
                %params.diameter = self.diameter;
                evalin('caller',['paramsStruct.intensity = ' num2str(self.intensity) ';']); %save the actual diameter used
            end
            
            %multiple integrated intensities
            if length(self.integratedIntensity) > 1
                if self.callCounter == 1 %if init
                    R = ceil(rand(1,nItems).*length(self.integratedIntensity));
                    integratedIntensityList = self.integratedIntensity(R);
                    self.integratedIntensityList = integratedIntensityList;
                end
                %use the next value in the list
                self.integratedIntensity = integratedIntensityList(self.callCounter);
                %params.diameter = self.diameter;
                evalin('caller',['paramsStruct.integratedIntensity = ' num2str(self.integratedIntensity) ';']); %save the actual diameter used
            end
            
            
            %make the texture
            radius = round(self.diameter/2);
            
            %check radius
            if self.Xpos-radius <= 0 || self.Xpos+radius > self.screenX || ...
                    self.Ypos-radius <= 0 || self.Ypos+radius > self.screenY
                disp('Spot Position error');
                self.initerror = 1;
                return;
            end
            
            %make the spot and save the texture
            M = ones(self.screenX,self.screenY).*self.spatial_meanLevel;
            square = ones(2*radius+1,2*radius+1).*self.spatial_meanLevel;
            X = -radius:radius;
            Y = -radius:radius;
            Xmat = repmat(X,2*radius+1,1);
            Ymat = repmat(Y,2*radius+1,1)';
            circlePoints = sqrt(Xmat.^2 + Ymat.^2)<radius;
            
            if ~self.integratedIntensity
                square(circlePoints) = self.intensity;
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
                square(circlePoints) = pixelValue;
            end
            
            M(self.Xpos-radius:self.Xpos+radius, self.Ypos-radius:self.Ypos+radius) = square;
            self.texturePtr = Screen('MakeTexture', self.windPtr, M');
            
            
        end
                
        function drawParams = nextStimFrame(self)
            drawParams.tex = self.texturePtr;
        end
        
    end
end
