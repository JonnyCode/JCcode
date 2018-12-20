classdef PairedSpotStimulus2 < SpatialStimSuper
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
        positionList
        intensityList
        integratedIntensityList
        pairSpots
        gammaTable        
        Xpos1
        Ypos1
        Xpos2
        Ypos2
    end
    
    methods
        
        function self = PairedSpotStimulus2(windPtr, params)
            self.initSuper(windPtr, params);
            
            persistent diameterList
            persistent intensityList
            persistent integratedIntensityList
            persistent positionList %X1, Y1, X2, Y2 (zeros for singles)

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
            
            %multiple X,Y positions
            if self.callCounter == 1 %if init
                if length(self.Xpos) ~= length(self.Ypos)
                    error('Xpos and Ypos must be the same length');
                end                                
                nSingles = length(self.Xpos)
                nPairs = nchoosek(nSingles,2);              
                pairIndList = combnk(1:nSingles,2)
                R = ceil(rand(1,nItems).*(nSingles+nPairs));
                positionList = zeros(nItems,4);
                for i=1:nItems
                   if R(i)<=nSingles
                      positionList(i,1) = params.Xpos(R(i));
                      positionList(i,2) = params.Ypos(R(i));
                      positionList(i,3) = 0;
                      positionList(i,4) = 0;
                   else
                      positionList(i,1) = params.Xpos(pairIndList(R(i)-nSingles,1));
                      positionList(i,2) = params.Ypos(pairIndList(R(i)-nSingles,1));
                      positionList(i,3) = params.Xpos(pairIndList(R(i)-nSingles,2));
                      positionList(i,4) = params.Ypos(pairIndList(R(i)-nSingles,2));
                   end
                end   
                self.positionList = positionList;
            end
            
            %self.positionList
            %use the next value in the list
            if positionList(self.callCounter,3) > 0 %if pair
                self.Xpos1 = positionList(self.callCounter,1);
                self.Ypos1 = positionList(self.callCounter,2);
                self.Xpos2 = positionList(self.callCounter,3);
                self.Ypos2 = positionList(self.callCounter,4);                
                self.pairSpots = 1;
                evalin('caller',['paramsStruct.Xpos1 = ' num2str(self.Xpos1) ';']); %save the actual pos used
                evalin('caller',['paramsStruct.Ypos1 = ' num2str(self.Ypos1) ';']);
                evalin('caller',['paramsStruct.Xpos2 = ' num2str(self.Xpos2) ';']);
                evalin('caller',['paramsStruct.Ypos2 = ' num2str(self.Ypos2) ';']);
                evalin('caller',['paramsStruct.pairSpots = ' num2str(1) ';']);
            else %if single
                self.Xpos1 = positionList(self.callCounter,1);
                self.Ypos1 = positionList(self.callCounter,2);
                self.pairSpots = 0;
                evalin('caller',['paramsStruct.Xpos1 = ' num2str(self.Xpos1) ';']); %save the actual pos used
                evalin('caller',['paramsStruct.Ypos1 = ' num2str(self.Ypos1) ';']);
                evalin('caller',['paramsStruct.pairSpots = ' num2str(0) ';']);
            end
            %make the texture
            radius = round(self.diameter/2);
            M = ones(self.screenX,self.screenY).*self.spatial_meanLevel;
            
            %keyboard;
            
            %self.pairSpots
            if self.pairSpots %pair of spots
                %check radius
                if self.Xpos1-radius <= 0 || self.Xpos1+radius > self.screenX || ...
                        self.Ypos1-radius <= 0 || self.Ypos1+radius > self.screenY || ...
                        self.Xpos2-radius <= 0 || self.Xpos2+radius > self.screenX || ...
                        self.Ypos2-radius <= 0 || self.Ypos2+radius > self.screenY
                    %error
                    disp('pair radius error')
                    self.initerror = 1;
                    return;
                end
                
                %make the spots and save the texture
                
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
                    disp(['Ivalue: ' num2str(Ivalue)]);
                    if Ivalue < 0 || Ivalue > 1
                        error(['Ivalue of ' num2str(Ivalue) 'is out of range']);
                    end
                    pixelValue = gammaCorrect(Ivalue,self.gammaTable);
                    disp(['PixelValue: ' num2str(pixelValue)]);
                    evalin('caller',['paramsStruct.intensity_raw=' num2str(pixelValue) ';']);
                    square(circlePoints) = pixelValue;
                end
                
                M(self.Xpos1-radius:self.Xpos1+radius, self.Ypos1-radius:self.Ypos1+radius) = square;
                M(self.Xpos2-radius:self.Xpos2+radius, self.Ypos2-radius:self.Ypos2+radius) = square;
                
            else %single spot
                
                %check radius
                if self.Xpos1-radius <= 0 || self.Xpos1+radius > self.screenX || ...
                        self.Ypos1-radius <= 0 || self.Ypos1+radius > self.screenY                    
                    disp('single radius error')
                    self.initerror = 1;
                    return;
                end
                
                %make the spot and save the texture
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
                    disp(['Ivalue: ' num2str(Ivalue)]);
                    if Ivalue < 0 || Ivalue > 1
                        error(['Ivalue of ' num2str(Ivalue) 'is out of range']);
                    end
                    pixelValue = gammaCorrect(Ivalue,self.gammaTable);
                    disp(['PixelValue: ' num2str(pixelValue)]);
                    evalin('caller',['paramsStruct.intensity_raw=' num2str(pixelValue) ';']);
                    square(circlePoints) = pixelValue;
                end
                
                M(self.Xpos1-radius:self.Xpos1+radius, self.Ypos1-radius:self.Ypos1+radius) = square;
                
            end %if pair
            
           % imagesc(M');
           % pause;
            self.texturePtr = Screen('MakeTexture', self.windPtr, M');
                        
        end
        
        
        function drawParams = nextStimFrame(self)
            drawParams.tex = self.texturePtr;
        end
        
    end
end
