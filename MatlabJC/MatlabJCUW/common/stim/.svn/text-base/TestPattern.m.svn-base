classdef TestPattern < SpatialStimSuper
    properties
        intensity
        diameter
        Xpos1
        Ypos1
        Xpos2
        Ypos2
        Xpos3
        Ypos3
    end
    properties (Hidden = true)
        gammaTable
    end

    methods

        function self = TestPattern(windPtr, params)
            self.initSuper(windPtr, params);

            self.Xpos1 = params.Xpos1;
            self.Ypos1 = params.Ypos1;
            self.Xpos2 = params.Xpos2;
            self.Ypos2 = params.Ypos2;
            self.Xpos3 = params.Xpos3;
            self.Ypos3 = params.Ypos3;
            self.diameter = params.diameter;
            self.gammaTable = params.gammaTable;
            self.intensity = gammaCorrect(params.intensity,self.gammaTable);
            evalin('caller',['paramsStruct.intensity_raw=[' num2str(self.intensity) '];']);

            %make the texture
            radius = round(self.diameter/2);
            M = ones(self.screenX,self.screenY).*self.spatial_meanLevel;

            %keyboard;

            %check radius
            if self.Xpos1-radius <= 0 || self.Xpos1+radius > self.screenX || ...
                    self.Ypos1-radius <= 0 || self.Ypos1+radius > self.screenY || ...
                    self.Xpos2-radius <= 0 || self.Xpos2+radius > self.screenX || ...
                    self.Ypos2-radius <= 0 || self.Ypos2+radius > self.screenY || ...
                    self.Xpos3-radius <= 0 || self.Xpos3+radius > self.screenX || ...
                    self.Ypos3-radius <= 0 || self.Ypos3+radius > self.screenY
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

            square(circlePoints) = self.intensity;

            M(self.Xpos1-radius:self.Xpos1+radius, self.Ypos1-radius:self.Ypos1+radius) = square;
            M(self.Xpos2-radius:self.Xpos2+radius, self.Ypos2-radius:self.Ypos2+radius) = square;
            M(self.Xpos3-radius:self.Xpos3+radius, self.Ypos3-radius:self.Ypos3+radius) = square;

            self.texturePtr = Screen('MakeTexture', self.windPtr, M');
        end %contructor


        function drawParams = nextStimFrame(self)
            drawParams.tex = self.texturePtr;
        end

    end
end
