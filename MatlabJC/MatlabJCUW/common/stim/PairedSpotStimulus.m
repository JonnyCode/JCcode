classdef PairedSpotStimulus < SpatialStimSuper
    properties
        Xpos
        Ypos
        diameter %can be a list of diameters to choose from
        intensity
        integratedIntensity %total intensity of spot, 0 to ignore and use intensity instead
        Npairs  %number of different pairs to sample
    end
    properties (Hidden = true)
        %persistent properties
        diameterList
        XposListSingles
        YposListSingles
        XposListPairs
        YposListPairs
        pairType
        pairSpots
        gammaTable
        Xpos1
        Ypos1
        Xpos2
        Ypos2
    end
    
    methods
        
        function self = PairedSpotStimulus(windPtr, params)
            self.initSuper(windPtr, params);
            
            persistent diameterList
            persistent XposListSingles
            persistent YposListSingles
            persistent XposListPairs
            persistent YposListPairs
            persistent pairType
            
            self.Xpos = params.Xpos;
            self.Ypos = params.Ypos;
            self.diameter = params.diameter;
            self.gammaTable = params.gammaTable;
            self.intensity = gammaCorrect(params.intensity,self.gammaTable);
            evalin('caller',['paramsStruct.intensity_raw=' num2str(self.intensity) ';']);
            self.integratedIntensity = params.integratedIntensity;
            self.Npairs = params.Npairs; %number of different pairs to sample
            Nsingles = length(self.Xpos);
            singlePairRatio = Nsingles/(Nsingles+self.Npairs);
            
            
            nItems = params.NumCycles; %may change name of this param
            
            
            
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
            if self.callCounter == 1 %if init
                if length(self.Xpos) ~= length(self.Ypos)
                    error('Xpos and Ypos must be the same length');
                end
                %get Npairs pairs first
                R1 = ceil(rand(1,self.Npairs*2).*length(self.Xpos));
                R2 = ceil(rand(1,self.Npairs*2).*length(self.Xpos));
                pairInd = find(R1~=R2);
                pairInd = pairInd(1:self.Npairs);
                Rpairs = [R1(pairInd); R2(pairInd)];
                
                pairType = (rand(1,nItems) > singlePairRatio);
                Rsingles = ceil(rand(1,Nsingles).*length(self.Xpos));
                
                XposListSingles = self.Xpos(Rsingles);
                YposListSingles = self.Ypos(Rsingles);
                XposListPairs = self.Xpos(Rpairs);
                YposListPairs = self.Ypos(Rpairs);
                self.XposListSingles = XposListSingles;
                self.YposListSingles = YposListSingles;
                self.XposListPairs = XposListPairs;
                self.YposListPairs = YposListPairs;
                self.pairType = pairType;
            end
            %use the next value in the list
            if pairType(self.callCounter) %if pair
                rand('seed',self.callCounter);
                R = ceil(rand*self.Npairs);
                self.Xpos1 = XposListPairs(1,R);
                params.Xpos1 = self.Xpos;
                self.Ypos1 = YposListPairs(1,R);
                params.Ypos1 = self.Ypos;
                
                self.Xpos2 = XposListPairs(2,R);
                params.Xpos2 = self.Xpos;
                self.Ypos2 = YposListPairs(2,R);
                params.Ypos2 = self.Ypos;
                self.pairSpots = 1;
                evalin('caller',['paramsStruct.Xpos1 = ' num2str(self.Xpos1) ';']); %save the actual pos used
                evalin('caller',['paramsStruct.Ypos1 = ' num2str(self.Ypos1) ';']);
                evalin('caller',['paramsStruct.Xpos2 = ' num2str(self.Xpos2) ';']);
                evalin('caller',['paramsStruct.Ypos2 = ' num2str(self.Ypos2) ';']);
                evalin('caller',['paramsStruct.pairSpots = ' num2str(1) ';']);
            else %if single
                rand('seed',self.callCounter);
                R = ceil(rand*Nsingles);
                self.Xpos = XposListSingles(R);
                params.Xpos = self.Xpos;
                self.Ypos = YposListSingles(R);
                params.Ypos = self.Ypos;
                self.pairSpots = 0;
                evalin('caller',['paramsStruct.Xpos = ' num2str(self.Xpos) ';']); %save the actual pos used
                evalin('caller',['paramsStruct.Ypos = ' num2str(self.Ypos) ';']);
                evalin('caller',['paramsStruct.pairSpots = ' num2str(0) ';']);
            end
            
            
            %multiple intensities
            
            %multiple integrated intensities
            
            %make the texture
            radius = round(self.diameter/2);
            M = ones(self.screenX,self.screenY).*self.spatial_meanLevel;
            
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
                if self.Xpos-radius <= 0 || self.Xpos+radius > self.screenX || ...
                        self.Ypos-radius <= 0 || self.Ypos+radius > self.screenY                    
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
                
                M(self.Xpos-radius:self.Xpos+radius, self.Ypos-radius:self.Ypos+radius) = square;
                
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
