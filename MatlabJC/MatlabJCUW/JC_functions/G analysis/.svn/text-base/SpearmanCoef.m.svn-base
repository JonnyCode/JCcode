function [Correlation,pvalue] = SpearmanCoef(signalX,signalY) ;

% this function calculates the Spearman Rank-Order Correlation Coefficient
% and indicates its significance 

% JC 1/10/08

% example signals
% X = normrnd(10,5,1,20000) ; % guassian
% Y = normrnd(10,5,1,20000) ; % guassian
% C = normrnd(0,0,1,20000) ; % guassian
% signalX = X+C; % matrix of signals where each row in siganlX should be correlated with a the same row in SignalY
% signalY = Y+C;

% rank both signals 
rankedX = ranker(signalX) ; 
rankedY = ranker(signalY) ;

% means of the ranked signals
meanRankedX = nanmean(rankedX) ; 
meanRankedY = nanmean(rankedY) ; 

% residual of the ranked signals
resdlRankedX = rankedX - repmat(meanRankedX,size(rankedX,1),1) ;
resdlRankedY = rankedY - repmat(meanRankedY,size(rankedY,1),1) ;

for a = 1:size(resdlRankedX,1) ; % for each signal to be correlated ...
    for timeshift = 1:length(resdlRankedX(a,:))+length(resdlRankedY(a,:))-1 ; % for each timeshift in the cross correlation ...
        if timeshift<=(length(resdlRankedX(a,:))+length(resdlRankedY(a,:)))/2 % if the timeshift is less than the middle shift than... 
            shiftedX = [NaNs(1,length(resdlRankedY(a,:))-timeshift),resdlRankedX(a,:)] ; % add NaNs so that the signals are equal in length 
            shiftedY = [resdlRankedY(a,:),NaNs(1,length(resdlRankedX(a,:))-timeshift)] ;
            numerator = nansum(shiftedX.*shiftedY) ;
            denominator = sqrt(nansum(shiftedX.^2))*sqrt(nansum(shiftedY.^2)) ; 
            Correlation(a,timeshift) = numerator/denominator ;
        else   
            shiftedX = [resdlRankedX(a,:),NaNs(1,timeshift-length(resdlRankedY(a,:)))] ; % add NaNs so that the signals are equal in length 
            shiftedY = [NaNs(1,timeshift-length(resdlRankedX(a,:))), resdlRankedY(a,:)] ;
            numerator = nansum(shiftedX.*shiftedY) ;
            denominator = sqrt(nansum(shiftedX.^2))*sqrt(nansum(shiftedY.^2)) ; 
            Correlation(a,timeshift) = numerator/denominator ;
        end
    end
    NoSlideCorr = Correlation(a,(length(resdlRankedX(a,:))+length(resdlRankedY(a,:)))/2) ; % the correlation without a time shift
    peakCorr = max(Correlation(a,:)) ; % the highest correlation between the two signals
    signifigance(a) = peakCorr*sqrt((length(resdlRankedX(a,:))-2)/(1-(peakCorr^2))) ; % the significance of the correlation = r*sqrt((N-2)/(1-r^2))    
    pvalue(a) = 2*(1-tcdf(signifigance(a),min(length(resdlRankedX),length(resdlRankedY)))) ; % probability that by chance, the absolute t value from the from student t-distribution is greater then the "signifigance" value calculated above     
end % NOTE - a pvalue below .05 would indicate signifigant correlation between the two signals

