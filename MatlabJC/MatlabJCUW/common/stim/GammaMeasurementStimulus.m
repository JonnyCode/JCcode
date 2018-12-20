classdef GammaMeasurementStimulus < SpatialStimSuper
    properties
        intensity_raw
        square_Width
    end
    
    methods
        function self = GammaMeasurementStimulus(windPtr, params)
            self.intensity_raw = params.intensity_raw;
            self.square_Width = params.square_Width;
                        
            self.initSuper(windPtr, params); 
            
            M=zeros(600,800) ;
            if self.intensity_raw==0 %for running a loop to go through all intensity values
                M(((600-self.square_Width)/2)+1:((600-self.square_Width)/2)+self.square_Width,((800-self.square_Width)/2)+1:((800-self.square_Width)/2)+self.square_Width) = (self.callCounter -1);         
            else
                M(((600-self.square_Width)/2)+1:((600-self.square_Width)/2)+self.square_Width,((800-self.square_Width)/2)+1:((800-self.square_Width)/2)+self.square_Width) = self.intensity_raw; 
            end
            self.texturePtr = Screen('MakeTexture', self.windPtr, M);   
        end
        
        function drawParams = nextStimFrame(self)
            drawParams.tex = self.texturePtr;            
        end
        
    end
    
end

