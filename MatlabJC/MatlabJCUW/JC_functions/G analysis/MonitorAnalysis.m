function [framerate_mean,framerate_range,startAmp] = MonitorAnalysis(Monitor, NumTrials, SI) 

% 1/11 JC adapted from GS
% 2/2/11 modified to deal with nans

for a = 1:NumTrials ;
    numPnts = find(isnan(Monitor(a,:))~=1,1,'last') ; % if nans are tacked to the end, ignore them
    
    MonitorLP = lowPassFilter(Monitor(a,1:numPnts),1/SI(1),100) ;
    MonitorLP = MonitorLP ;
    MonitorLPNorm = MonitorLP/max(MonitorLP(:)) ;
    MonitorLPNormShift = circshift(MonitorLPNorm,[0,1]) ;

    time = [SI:SI:numPnts*SI] ;
    
    CrossPnts = find(MonitorLPNorm>.5 & MonitorLPNormShift<=.5)  ;
    diffCrossPnts = diff(CrossPnts) ;
    
    timePnts = time(CrossPnts(2:end)) ;
    framerate =  60./diff(timePnts) ; % 60 frames per cross thus framerate is in frames/sec
    framerate_mean(a) = mean(framerate) ;
    framerate_range(a) = max(framerate)-min(framerate) ;
end

startAmp = MonitorLPNorm(:,1) ;