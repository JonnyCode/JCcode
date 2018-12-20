classdef ThreeDdotStimulus < SpatialStimSuper
    properties %separate into private and public (params needed from IGOR)
        maskRadius      % radius of mask (pix)
        xPosMean           % mean of x vector (position in respect to screen 0=left & 1=far right)
        yPosMean           % mean of y vector (position in pix)
        zPosMean           % mean of z vector (size in pix to simulate 3rd dimension)
        xPosStd            % std of x vector (drift of position in pix)
        yPosStd
        zPosStd            % std of z vector (drift of size in pix to simulate drift in 3rd dimension)
        xPosFreqCut        % highest freqency of x vector (pix/frame)
        yPosFreqCut        
        zPosFreqCut
        NoDrift         % 0=normal drift, 1=spot moves along across screen without changing direction
        
    end
    properties (Hidden = true)
        randStream
        gammaTable
        BarPts
        BarTexture
        BlankTexture
        maskLength
    end
    methods
        
        function self = ThreeDdotStimulus(windPtr, params)
            self.initSuper(windPtr, params) ;                
            
            maskRadius = params.maskRadius ;
            xPosMean = params.xPosMean ;
            yPosMean = params.yPosMean ;
            zPosMean = params.zPosMean ;
            xPosStd = params.xPosStd ;
            yPosStd = params.yPosStd ;
            zPosStd = params.zPosStd ;
            xPosFreqCut = params.xPosFreqCut ;
            yPosFreqCut = params.yPosFreqCut ;
            zPosFreqCut = params.zPosFreqCut ;
            NoDrift = params.NoDrift ;
            
            % make each position/size vector
            FakeSampleRate = 10000 ;
            
            x = normrnd(xPosMean,xPosStd,1,self.spatial_stimpts) ;
            y = normrnd(yPosMean,yPosStd,1,self.spatial_stimpts) ;
            z = normrnd(zPosMean,zPosStd,1,self.spatial_stimpts) ;
            
            lpFilteredSignal = lowPassFilter(signal, samplerate, freqcutoff) ;
            x = lowPassFilter(x, FakeSampleRate, freqcutoff) ;
            
            
            % FUNCTION INCOMPLETE
            
            
            
            
            
            
            