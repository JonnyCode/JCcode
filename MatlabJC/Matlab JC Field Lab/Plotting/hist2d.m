function hist2d = hist2d(vector1,vector2,vector1Bins,vector2Bins) 

% vectorBins are center bins

hist2d = zeros(length(vector1Bins),length(vector2Bins)) ; % prep hist

for a=1:length(vector1) ; % for each point in 1
    [m,r] = min((vector1(a)-vector1Bins).^2) ;
    [m,c] = min((vector2(a)-vector2Bins).^2) ;
    
    hist2d(r,c) = hist2d(r,c)+1 ;
end
        
        
        
        