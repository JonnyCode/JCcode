classdef  TextureStimulus < SpatialStimSuper
    properties
        size %pixels (edge length of square)
        corrFactor %s.d. of gaussian to convolve (pixles, can be non-integral)      
        smoothWin
        edgeSmoothFactor
        contrast %0 to 1
    end
    properties (Hidden = true)
        %persistent properties
        sizeList
        corrFactorList   
        contrastList
        gammaTable
    end
    
    methods
        
        function self = TextureStimulus(windPtr, params)
            self.initSuper(windPtr, params);
                       
            persistent sizeList
            persistent corrFactorList
            persistent contrastList
           
            self.smoothWin = params.smoothWin;
            self.edgeSmoothFactor = params.edgeSmoothFactor;
            self.size = params.size;
            self.corrFactor = params.corrFactor;
            self.contrast = params.contrast;
            self.gammaTable = params.gammaTable;   
            
            nItems = params.NumCycles;
            
            %multiple sizes
            if length(self.size) > 1
                if self.callCounter == 1 %if init
                    R = ceil(rand(1,nItems).*length(self.size));
                    sizeList = self.size(R);
                    self.sizeList = sizeList;
                end
                %use the next value in the list
                self.size = sizeList(self.callCounter);
                evalin('caller',['paramsStruct.size = ' num2str(self.size) ';']); %save the actual size used
            end
            
            %multiple corrFactors 
            if length(self.corrFactor) > 1
                if self.callCounter == 1 %if init
                    R = ceil(rand(1,nItems).*length(self.corrFactor));
                    corrFactorList = self.corrFactor(R);
                    self.corrFactorList = corrFactorList;
                end
                %use the next value in the list
                self.corrFactor = corrFactorList(self.callCounter);
                evalin('caller',['paramsStruct.corrFactor = ' num2str(self.corrFactor) ';']); %save the actual size used
            end
            
            %multiple contrasts 
            if length(self.contrast) > 1
                if self.callCounter == 1 %if init
                    R = ceil(rand(1,nItems).*length(self.contrast));
                    contrastList = self.contrast(R);
                    self.contrastList = contrastList;
                end
                %use the next value in the list
                self.contrast = contrastList(self.callCounter);
                evalin('caller',['paramsStruct.contrast = ' num2str(self.contrast) ';']); %save the actual size used
            end
            
            M = ones(params.screenX,params.screenY)*self.spatial_meanLevel;
            
            startX = params.screenX/2-round(self.size/2);
            startY = params.screenY/2-round(self.size/2);
            %generate the texture
            T = generateTexture(self.size,self.corrFactor,self.contrast);
            
            %gamma correct each row
            for i=1:self.size
               T(i,:) = gammaCorrect(T(i,:), self.gammaTable);
            end
            
            %put texture in window
            M(startX:startX+self.size-1, startY:startY+self.size-1) = T;
                      
%             %smooth edges
%             winL = self.smoothWin;    
%             sFact = self.edgeSmoothFactor;
%             W = linspace(0,1,winL);
% %            W = gausswin(winL)';
%             Winv = 1-W;
%             F1 = W./sum(W);
%             F2 = Winv./sum(Winv);
%             
%             Gwin = 8;
%             win = fspecial('gaussian',Gwin,sFact);
% 
%             %LR
% %             leftEdgeInd = [startY-winL:startY+winL];
% %             rightEdgeInd = [startY+self.size-winL:startY+self.size+winL];            
% %             for i=1:size(M,1)
% %                 if i>=startX && i<startX+self.size
% %                     CpartL = conv2(M(i,leftEdgeInd),F1,'valid');
% %                     CpartR = conv2(M(i,rightEdgeInd),F2,'valid');
% %                     M(i,startY-winL/2-1:startY+winL/2) = CpartL;
% %                     M(i,startY+self.size-winL/2-1:startY+self.size+winL/2) = CpartR;
% %                 end
% %             end        
% 
%             %LR
%             leftEdgeInd = [startY-winL:startY+winL];
%             rightEdgeInd = [startY+self.size-winL:startY+self.size+winL];
%             CpartL = conv2(M(:,leftEdgeInd),win,'valid');
%             CpartR = conv2(M(:,rightEdgeInd),win,'valid');
%             %size(CpartL)
%             %size(M(Gwin/2:end-Gwin/2,startY-winL+Gwin/2-1:startY+winL-Gwin/2))
%             M(Gwin/2:end-Gwin/2,startY-winL+Gwin/2-1:startY+winL-Gwin/2) = CpartL;
%             M(Gwin/2:end-Gwin/2,startY-winL+Gwin/2-1+self.size:startY+winL-Gwin/2+self.size) = CpartR;
% 
%             %UD
%             topEdgeInd = [startX-winL:startX+winL];
%             bottomEdgeInd = [startX+self.size-winL:startX+self.size+winL];
%             CpartT = conv2(M(topEdgeInd,:),win,'valid');
%             CpartB = conv2(M(bottomEdgeInd,:),win,'valid');
%             M(startX-winL+Gwin/2-1:startX+winL-Gwin/2,Gwin/2:end-Gwin/2) = CpartT;
%             M(startX-winL+Gwin/2-1+self.size:startX+winL-Gwin/2+self.size,Gwin/2:end-Gwin/2) = CpartB;
% 
%             
% %             %UD
% %             topEdgeInd = [startX-winL:startX+winL];
% %             bottomEdgeInd = [startX+self.size-winL:startX+self.size+winL];            
% %             for i=1:size(M,2)
% %                 if i>=startY && i<startY+self.size
% %                     CpartT = conv2(M(topEdgeInd,i),F1','valid');
% %                     CpartB = conv2(M(bottomEdgeInd,i),F2','valid');
% %                     M(startX-winL/2-1:startX+winL/2,i) = CpartT;
% %                     M(startX+self.size-winL/2-1:startX+self.size+winL/2,i) = CpartB;
% %                 end
% %             end
            
            %make openGL texture
            self.texturePtr = Screen('MakeTexture', self.windPtr, M');
                        
        end
                
        function drawParams = nextStimFrame(self)
            drawParams.tex = self.texturePtr;
        end
        
    end
end