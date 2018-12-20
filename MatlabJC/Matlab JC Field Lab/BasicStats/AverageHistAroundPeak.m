function histFolded = AverageHistAroundPeak(histVect,varargin)

% this function will average a vector from values around its peak assuming a 360 degree x-axis.
% JC 2018-01-17

if length(varargin)==0 ; % if not given a point to average around
    [m,mi] = max(histVect) ; % find first max value
else
    mi = varargin{1} ; % use given index as center
end

histVectShift = circshift(histVect,[0,-mi+1]) ; % circ shift so max is first point
histFold = (histVectShift(2:end)+ fliplr(histVectShift(2:end)))/2 ; % average across points from peak
histFolded = [m,histFold(1:ceil(length(histFold)/2))] ; % average of hist values surrounding peak (peak is first point)
