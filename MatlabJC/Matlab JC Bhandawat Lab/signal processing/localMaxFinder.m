function [localMaxPnts,localMaxValues] = localMaxFinder(signal) ;

% this function will find the local maximum points of a vector "signal"

signal_diff = diff(signal) ;
localMaxPnts = nan(1,round(length(signal_diff)/2)) ;
localMaxValues = nan(1,round(length(signal_diff)/2)) ;

r=0 ;
for a=2:length(signal_diff) ;
    if signal_diff(a)<=0 && signal_diff(a-1)>0 ;
        r=r+1 ;
        localMaxPnts(r) = a ;
        localMaxValues(r) = signal(a) ;
    end
end

if r>0 ;
    localMaxPnts = localMaxPnts(1:r) ;
    localMaxValues = localMaxValues(1:r) ;
else
    localMaxPnts = nan ;
    localMaxValues = nan ;
end