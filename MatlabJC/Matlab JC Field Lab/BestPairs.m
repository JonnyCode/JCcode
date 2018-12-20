function Pairs = BestPairs(CorrMat,Threshold) ;

% function to create best X,Y pairs whose correlations exceed a threshold

% "CorrMat" length(row)>=length(columns)

% JC 6/11/15

HighN = size(CorrMat,1) ;
M = nans(HighN,1) ;
modMat = CorrMat.*(modMat>Threshold) ; % find all possible matches

while sum(isnan(M)>0) ; % while there is still a possible match
    for i=1:HighN ; % for each possible
        if sum(modMat(i,:))>0 ; % if there are possible matches
            [v,j] = max(modMat(i,:)) ; % find the best possible match 
            if sum(M==j)<1 ; % if this cell has not been used as a match
                M(i) = j ; % that is the match
            else % if this cell has been used as a match
                if modMat(i,j)<modMat(M==j,j) ; % if the old match is better
                    modMat(i,j) = 0 ; % not possibility anymore
                else % if the new match is better
                    M(M==j) = nan ; % keep searching for old match
                    modMat(M==j,j)= 0 ; % that one is no longer possible for old match 
                    M(i) = j ; % give it to the new match
                end
            end
        else
            M(i) = [] ; % no match possible
        end
    end
end
                        
