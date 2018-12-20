function Output = NonLinFilter2(InputBins1,InputBins2,Input1,Input2,Nonlinearity) ;

% this function will implement the 3d nonlinearity generated in
% NonLinFilterFinder2.m on a given input.  This function will not
% interpolate outside of its range.

% JC 8/14/10

for a = 1:numel(Input1) ;
    if Input1(a)>=min(InputBins1) & Input1(a)<=max(InputBins1) & Input2(a)>=min(InputBins2) & Input2(a)<=max(InputBins2) ; % if your within range of bins
        Output(a) = interp2(InputBins1,InputBins2,Nonlinearity,Input1(a),Input2(a)) ; % assume linearity between bins
    elseif Input1(a)>=min(InputBins1) & Input1(a)<=max(InputBins1) ;                                                    % if your within range of Input bins 1 only
        InputBin1_nearest = InputBins1(find(min(InputBins1-Input1(a)))) ; % the nearest relavent bin within range
        Output(a) = interp(InputBins2,Nonlinearity(InputBin1_nearest,:),Input2(a),'linear','extrap') ; % linear extrapolation down row
    elseif Input2(a)>=min(InputBins2) & Input2(a)<=max(InputBins2) ;                                                    % if your within range of Input bins 2 only
        InputBin2_nearest = InputBins2(find(min(InputBins2-Input2(a)))) ; % the nearest relavent bin within range
        Output(a) = interp(InputBins1,Nonlinearity(InputBin2_nearest,:),Input1(a),'linear','extrap') ; % linear extrapolation down column
    else                                                                                                                % if your outside the range of both
        Output(a) = nan ; % THIS SHOULD BE FIXED SO THAT FUNCTION INTERPOLATES THESE VALUES
    end
            
        
end
    

