function params = tc_params_finder(tc,tc_time) 

% tc = time course of STA matrix with rows as idividual tcs (spike time is first point)
% tc_time = time vector of STA time course (looking for [0:-t])

% JC 5/20/15

for a= 1:size(tc,1) ; % for each tc
    tca = tc(a,:) ;
    [tcmin,mini] = min(tca) ;
    [tcmax,maxi] = max(tca) ;
    
    if tcmin<0 && tcmax>0 ; % if there are both positive and negative values
        firstPeaki = min(mini,maxi) ;
        secondPeaki = max(mini,maxi) ;
        
        params.firstPeakt(a) = tc_time(firstPeaki) ; 
        params.secondPeakt(a) = tc_time(secondPeaki) ; 
        params.firstPeak(a) = tca(firstPeaki) ;
        params.secondPeak(a) = tca(secondPeaki) ;
        params.PeakRatio(a) = tca(firstPeaki)/tca(secondPeaki) ;
        
        params.zeroCrosst(a) = interp1(tca(firstPeaki:secondPeaki),tc_time(firstPeaki:secondPeaki),0,'linear') ;
        
        zeroCrossi = find(tc_time>params.zeroCrosst(a),1,'last') ;
        params.firstPeakArea(a) = trapz(tca(1:zeroCrossi))*diff(tc_time(1:2)) ; 
        params.secondPeakArea(a) = trapz(tca(zeroCrossi+1:end))*diff(tc_time(1:2)) ; 
        params.DoT(a) = abs(params.firstPeakArea(a)+params.secondPeakArea(a))/(abs(params.firstPeakArea(a))+abs(params.secondPeakArea(a))) ;  
        
    else
        params.firstPeakt(a) = nan ; 
        params.secondPeakt(a) = nan ; 
        params.firstPeak(a) = nan ;
        params.secondPeak(a) = nan ;
        params.PeakRatio(a) = nan ;    
        params.zeroCrosst(a) = nan ;
        params.firstPeakArea(a) = nan ; 
        params.secondPeakArea(a) = nan ; 
        params.DoT(a) = nan ; 
    end
      
end
        
