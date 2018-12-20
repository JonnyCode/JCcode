classdef SpatialGradingStimulus < SpatialStimSuper
    properties %separate into private and public (params needed from IGOR)
        BarWidth        % pixel width of half cycle of square wave
        Phase           % fraction of cycle to move square wave over
        spotRadius      % radius of mask
        amplitude       % change from mean
        CycleFrames     % number of frames during a cycle
        squareWaveSpace % 1 = square, else = sinewave
        squareWaveTime  % 1 = square, else = sinewave
        recWaves        % 0 = normal sine, 1 = only above wave, -1=only below waves
        sumWaves        % 1 = sum all waves and run through all + sum, 0 = don't
    end
    properties (Hidden = true)
        texture
        index

    end
    methods
        
        function self = SpatialGradingStimulus(windPtr, params)
            self.initSuper(windPtr, params) ;
            
            BarWidth = params.BarWidth ;
            Phase = params.Phase ;
            spotRadius = params.spotRadius ;  
            amplitude = params.amplitude ;
            self.CycleFrames = params.CycleFrames ;    
            squareWaveTime = params.squareWaveTime ;
            squareWaveSpace = params.squareWaveSpace ;
            recWaves = params.recWaves ;
            sumWaves = params.sumWaves ;
            
            % appropriate trial (if multiple inputs)
            self.index = rem(self.callCounter-1,length(BarWidth)+sumWaves)+1 ;
            
            if self.rand_behavior==1 ;
                R = randperm(length(BarWidth)+sumWaves);
                self.index = R(self.index) ;
            end
            
            
            if self.index<=length(BarWidth) ; % if this is not a sum wave trial
                % make sine waves
                sineWave_space = sin((pi/BarWidth(self.index))*[1:self.screenX]+Phase(self.index)*pi/180) ;
                sineWave_time = sin((pi/(self.CycleFrames(self.index)/2)*[1:self.CycleFrames(self.index)])) ;

                % change into square waves if desired
                if squareWaveSpace(self.index) == 1 ; 
                    sineWave_space(sineWave_space>=0)=1 ;
                    sineWave_space(sineWave_space<0)=-1 ;
                end

                if squareWaveTime(self.index) == 1 ;
                    sineWave_time(sineWave_time>=0)=1 ;
                    sineWave_time(sineWave_time<0)=-1 ;
                end
                
                % mask
                if spotRadius(self.index)~=0 ;              
                    xmask = repmat([1:self.screenX],self.screenY,1)-(self.screenX/2) ;
                    ymask = repmat([1:self.screenY]',1,self.screenX)-(self.screenY/2) ;
                    maskPnts = sqrt(xmask.^2 + ymask.^2) > spotRadius(self.index) ;
                else 
                    maskPnts = [] ;
                end
                
                % get correct intensity
                sineWave_space = sineWave_space*amplitude(self.index) ; 
%keyboard
                % get one wave for every frame of cycle (in rows)
                Waves = repmat(sineWave_space,self.CycleFrames(self.index),1).*repmat(sineWave_time',1,self.screenX) + params.spatial_meanLevel ; % each spatial wave is mulitplied by temporal wave going from 1 to -1  
                
                % recitify waves if desired
                if recWaves(self.index) == 1 ;
                    Waves(Waves<params.spatial_meanLevel)=params.spatial_meanLevel ;
                elseif recWaves(self.index) == -1 ;
                    Waves(Waves>params.spatial_meanLevel)=params.spatial_meanLevel ;
                end
                
                % gamma correct
                Waves = gammaCorrect(Waves,params.gammaTable) ;

                % textures
                for a=1:self.CycleFrames(self.index) ;
                    textureMatrix = repmat(Waves(a,:),self.screenY,1) ;             
                    textureMatrix(maskPnts) =  self.spatial_meanLevel ;
                    self.texture{a} = Screen('MakeTexture', self.windPtr, textureMatrix) ;
                end
                
%                 % record proper parameters
                if self.index<=length(BarWidth) ; % if this is not a sum wave trial              
                    evalin('caller',['paramsStruct.BarWidth = ' num2str(BarWidth(self.index)) ';']) ;
                    evalin('caller',['paramsStruct.Phase = ' num2str(Phase(self.index)) ';']) ;
                    evalin('caller',['paramsStruct.spotRadius = ' num2str(spotRadius(self.index)) ';']) ;
                    evalin('caller',['paramsStruct.amplitude = ' num2str(amplitude(self.index)) ';']) ;
                    evalin('caller',['paramsStruct.CycleFrames = ' num2str(self.CycleFrames(self.index)) ';']) ;
                    evalin('caller',['paramsStruct.squareWaveSpace = ' num2str(squareWaveSpace(self.index)) ';']) ;
                    evalin('caller',['paramsStruct.squareWaveTime = ' num2str(squareWaveTime(self.index)) ';']) ;
                end

                
                
            else  % make sum of sinewaves
                Waves = zeros(max(self.CycleFrames),self.screenX) ;
                for a = 1:length(BarWidth) ; % make each sinewave
                    % make sine waves
                    sineWave_space = sin((pi/BarWidth(a))*[1:self.screenX]+Phase(a)*pi/180) ;
                    sineWave_time = sin((pi/(self.CycleFrames(a)/2)*[1:max(self.CycleFrames)])) ;

                    % change into square waves if desired
                    if squareWaveSpace(a) == 1 ; 
                        sineWave_space(sineWave_space>=0)=1 ;
                        sineWave_space(sineWave_space<0)=-1 ;
                    end

                    if squareWaveTime(a) == 1 ; %#ok<*PROP>
                        sineWave_time(sineWave_time>=0)=1 ;
                        sineWave_time(sineWave_time<0)=-1 ;
                    end
                    
                    % get correct intensity
                    sineWave_space = sineWave_space*amplitude(a) ; 
                
                    % get one wave for every frame of cycle (in rows) and add to previous waves
                    temp = repmat(sineWave_space,max(self.CycleFrames),1).*repmat(sineWave_time',1,self.screenX) ; % each spatial wave is mulitplied by temporal wave going from 1 to -1 
                    
                    % recitify waves if desired
                    if recWaves(a) == 1 ;
                        temp(temp<params.spatial_meanLevel)=params.spatial_meanLevel ;
                    elseif recWaves(a) == -1 ;
                        temp(temp>params.spatial_meanLevel)=params.spatial_meanLevel ;
                    end
                    
                    Waves = Waves + temp ;
                
                end
  %keyboard               
                Waves = gammaCorrect(Waves+params.spatial_meanLevel,params.gammaTable) ; 
            
                % mask
                if min(spotRadius)~=0 ;              
                    xmask = repmat([1:self.screenX],self.screenY,1)-(self.screenX/2) ;
                    ymask = repmat([1:self.screenY]',1,self.screenX)-(self.screenY/2) ;
                    maskPnts = sqrt(xmask.^2 + ymask.^2) > mean(spotRadius) ;
                else 
                    maskPnts = [] ;
                end

                % textures
                for a=1:max(self.CycleFrames) ;
                    textureMatrix = repmat(Waves(a,:),self.screenY,1) ;             
                    textureMatrix(maskPnts) =  self.spatial_meanLevel ;
                    self.texture{a} = Screen('MakeTexture', self.windPtr, textureMatrix) ;
                end
                %keyboard
            end            
        end
       
        
        function drawParams = nextStimFrame(self)
            frame = self.frameCounter-self.spatial_prepts ;
            
            if self.index<=length(self.CycleFrames) ; % if this is not a sum wave trial
                framePull = rem(frame-1,self.CycleFrames(self.index))+1 ;
            else
                framePull = rem(frame-1,max(self.CycleFrames))+1 ;
            end
                % pause
                drawParams.tex = self.texture{framePull};

        end
    end
end



