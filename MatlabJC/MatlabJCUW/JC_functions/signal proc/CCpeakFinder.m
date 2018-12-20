function CCpeak = CCpeakFinder(cc) ;

% this will find the peak deviation of a cross correlation or anything else
% that can have a positive of negative peak

CCpeak = cc(find(abs(cc)==max(abs(cc)),1)) ;

if isempty(CCpeak) ; % if cc is nans then peak should be nan
    CCpeak = nan ;
end

end
