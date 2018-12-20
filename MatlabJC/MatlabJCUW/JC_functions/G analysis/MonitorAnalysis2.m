function FrameStartPnts = MonitorAnalysis2(Monitor, NumTrials, SI) 

% 8/21/11 adapted from earlier version ('MonitorAnalysis.m') to get more
% accurate measure of frame start time.

UpSweepTime = .004 ; %sec :  if Monitor trace increases for this much time you are likely during the first frame 
UpSweepPnts = floor(UpSweepTime/SI(1)) ;

InterUpSweepTime_estimatedMin = .9 ; %sec ; likely minumum time between upsweeps (which should happen every 60 frames) 
InterUpSweepPnts_estimatedMin = floor(InterUpSweepTime_estimatedMin/SI(1)) ;

FrameStartPnts = cell(1,NumTrials) ;
for a = 1:NumTrials ;
    numPnts = find(isnan(Monitor(a,:))~=1,1,'last') ; % if nans are tacked to the end, ignore them
    MonitorDiff = diff(Monitor(a,1:numPnts)) ;
    
    FrameStartPnts{a} = nan(1,ceil(numPnts/InterUpSweepPnts_estimatedMin)) ;
    FrameStartPnts{a}(1) = 1 ; % The first point is attomatically considered the first upsweep point!
    MonitorIndex = 1 + InterUpSweepPnts_estimatedMin ;
    rnd = 2 ;
    while MonitorIndex<numPnts - UpSweepPnts ;
        if min(MonitorDiff(MonitorIndex:MonitorIndex+UpSweepPnts))>0 ; % if every point during the possible upsweep time is increasing than you are at the begining of an upsweep
            FrameStartPnts{a}(rnd) = MonitorIndex ;
            rnd = rnd + 1 ;
            MonitorIndex = MonitorIndex + InterUpSweepPnts_estimatedMin ; % jump about an upsweep forward
        else
            MonitorIndex = MonitorIndex + 1 ;
        end
    end
    
    FrameStartPnts{a} = floor(FrameStartPnts{a}(1:rnd-1)) ; % get rid of extras
    
end
