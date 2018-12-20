classdef AnalysisSuper < handle

    properties (Hidden = true)
        DataFilePath
        SpikeTrain
        fp        
        EpochData
        SampleInterval
        Time
    end
    methods
        function initAnalysisSuper(self, params)

            self.DataFilePath = params.dataFileName;
            % self.DataFilePath = ['/Volumes/' regexprep(params.dataFileName, ':', '/')]; %fix IGOR path issues
            %self.DataFilePath = '/Network/Servers/rieke-server.physiol.washington.edu/Volumes/Data/Users/Greg/acquisition/073109Bc1';

            
            [self.fp, error] = ITCInitializeAnalysis(1000000, self.DataFilePath); % get data file pointer            
        end
        
        function EpochDataPuller(self,epochNums,channel)  
            % epochNums = vector, Channel = 0 or 1            
            for a = 1:length(epochNums) ;
                [self.EpochData(a,:), error] = ITCReadEpoch(epochNums(a), channel, self.fp); % data into matrix channel 1
                
                [si, error] = ITCGetSamplingInterval(epochNums(a), self.fp); % sampling interval into vector
                self.SampleInterval(a) = si * 10^-6 ; % sample interval in sec
            end
        end
        
        
        function SpikeDetect(self,Data,threshold) % spike detection (threshold can be - or +) 
            
            DataShift = circshift(Data,[0,1]) ; % shift 
            
            self.SpikeTrain = zeros(size(Data)) ; % prep spike trains
            
            if threshold>0 ;
                self.SpikeTrain(Data>=threshold & DataShift<threshold) = 1 ; % find indicies where change in current is above threshold
            else
                self.SpikeTrain(Data<=threshold & DataShift>threshold) = 1 ;
            end
            
            self.SpikeTrain(:,1) = 0 ; % cannot dectect first point as spike
                       
        end
        
        function MakeTime(self,epochNum) % returns a vector of timepoints (all epochs better have the same length)
            [prep,error] = ITCGetStmPrePts(epochNum,0,0,self.fp);
            [stmp,error] = ITCGetStmPts(epochNum,0,0,self.fp);
            [taip,error] = ITCGetStmTailPts(epochNum,0,0,self.fp);
            self.Time = ((1:(prep+stmp+taip)) - prep - 1)*1e-3;  % hack
            
        end
        
    end
end
        