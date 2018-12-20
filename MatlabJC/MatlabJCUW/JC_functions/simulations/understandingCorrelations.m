T=5000 ;

for a=1:50 ;
    v1(a,:) = normrnd(normrnd(10,5,1,1),3,1,T) ; % random vectors with different means but same variance
end
v2 = circshift(v1,[0,300]) ; % correlated vector shifted in time

for a=1:50 ; % residuals
    v1r(a,:) = v1(a,:)-mean(v1) ;
    v2r(a,:) = v2(a,:)-mean(v2) ;
end


for a=1:50 ; % cross covariance
    cc(a,:) = xcorr(v1r(a,:),v2r(a,:),'biased') ; % mean of cross correlation for time dependant residuals
    ac1(a,:) = xcorr(v1r(a,:),'biased') ;
    ac2(a,:) = xcorr(v2r(a,:),'biased') ;
end

% averaged across time
MeanCov = cov(mean(v1,2),mean(v2,2)) 

cc_mean_mean = 2*mean(mean(cc)) 

% summed accross time
SumCov = cov(sum(v1,2),sum(v2,2)) 

cc_mean_sum = sum(mean(cc))*T ;

ac1_mean_sum = sum(mean(ac1))*T ;
ac2_mean_sum = sum(mean(ac2))*T ;

% corr coef
CorrCoef = corrcoef(sum(v1,2),sum(v2,2)) ;

ccCoef = cc_mean_sum/sqrt(ac1_mean_sum*ac2_mean_sum) ;







