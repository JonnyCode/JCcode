classdef StimulusStabilityAnalysis < AnalysisSuper
    properties
        epochNums        
    end
    
    methods
        
        function self = StimulusStabilityAnalysis(params)
            self.initAnalysisSuper(params) ;
            
            self.epochNums = str2num(params.epochNums) ;
            
            EpochDataPuller(self,self.epochNums,1) ;
            
            time = [.0001:.0001:length(self.EpochData)*.0001] ;
            
            Data = self.EpochData./max(self.EpochData(1:end)) ;
            shiftData = circshift(Data,[0,1]) ;
            
            for a=1:size(Data,1) ;

                i = find(Data(a,:)>=max(Data(1:end))/2 & shiftData(a,:)<max(Data(1:end))/2,1,'last') ;
                trig(a) = time(i) ;

            end
            
            meanTrig = mean(trig) 
            stdTrig = std(trig)   
            
            figure
            subplot(2,1,1)
            plot(time,Data)
            xlabel('time (sec)')
            ylabel('lightmeter normalized')
            
            subplot(2,1,2)
            plot(time,Data)
            hold on
            plot(trig,max(Data(1:end))/2,'ro')
            xlabel('time (sec)')
            ylabel('lightmeter normalized')
            axis([min(trig)-.25, max(trig)+.25,-1,2]) 
                     
            drawnow
            pause
            
        end
    end
end
