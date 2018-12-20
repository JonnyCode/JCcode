% 10/15/10 simulation to look at drift correcting adapting cell pairs

% parameters

for sample = 1:1000 ;

    numTrials = 100 ;
    adapt_Tau1 = 20 ; % adaptation tau
    adapt_Tau2 = 20 ; % adaptation tau

    corr = 0 ; 
    max_response = 100 ; % maximum average response
    min_response = 0 ; % minimim average response
    
    trials = [1:numTrials] ;

    % make signals in 2 cells
    signal1 = (max_response-min_response)*exp((-trials+1)/adapt_Tau1)+min_response ;
    signal2 = (max_response-min_response)*exp((-trials+1)/adapt_Tau2)+min_response ;

    % make noise

    noiseCommon = normrnd(0,sqrt(corr),1,numTrials) ; % common noise
    noiseInd1 = normrnd(0,sqrt(1-corr),1,numTrials) ;
    noiseInd2 = normrnd(0,sqrt(1-corr),1,numTrials) ;

    noise1 =(noiseInd1+noiseCommon).*sqrt(signal1) ; % scale noise such that variance=mean
    noise2 =(noiseInd2+noiseCommon).*sqrt(signal2) ;

    % noise1 =(noiseInd1+noiseCommon).*2 ; % scale noise such that variance is constant
    % noise2 =(noiseInd2+noiseCommon).*2 ;


    % add signal and noise
    measured1 = signal1 + noise1 ;
    measured2 = signal2 + noise2 ;


    % calculate noise correlation

    for a=2:numTrials-1 ; 
        res1(a-1) = measured1(a) - (measured1(a-1)+measured1(a+1))/2 ;
        res2(a-1) = measured2(a) - (measured2(a-1)+measured2(a+1))/2 ;
    end


%     for a=1:numTrials ; % standard mean calculation
%         res1(a) = measured1(a) - mean(measured1) ;
%         res2(a) = measured2(a) - mean(measured2) ;
%     end

    % correlation coef
    [cc,p] = corrcoef(res1',res2') ; 

    corrMeasured(sample) = cc(2,1) ; % no p value
   
    if p(2,1)<.005 ;
        corrMeasuredsig(sample) = cc(2,1) ;
    else
        corrMeasuredsig(sample) = 0 ;
    end
    
   % correlation coef shuffle

   res1_shuffle = nans(size(res1)) ;
   res1_shuffle(1:2:end) = res1(2:2:end) ;
   res1_shuffle(2:2:end) = res1(1:2:end) ;
    
    % correlation coef
    [cc,p] = corrcoef(res1_shuffle',res2') ; 

    corrMeasured_shuff(sample) = cc(2,1) ; % no p value

    if p(2,1)<.005 ;
        corrMeasuredsig_shuff(sample) = cc(2,1) ;
    else
        corrMeasuredsig_shuff(sample) = 0 ;
    end

end

corrMeasured_hist = hist(corrMeasured,[-1:.1:1]) ;
corrMeasuredsig_hist = hist(corrMeasuredsig,[-1:.1:1]) ;

corrMeasured_shuff_hist = hist(corrMeasured_shuff,[-1:.1:1]) ;
corrMeasuredsig_shuff_hist = hist(corrMeasuredsig_shuff,[-1:.1:1]) ;

figure
plot([-1:.1:1],corrMeasured_hist,'g')
hold on
plot([-1:.1:1],corrMeasuredsig_hist,'k--')
plot([-1:.1:1],corrMeasuredsig_shuff_hist,'r--')
xlabel('cc')
ylabel('number of observations')
legend('cc','cc sig','cc sig shuff')

figure
plot(corrMeasuredsig,corrMeasuredsig_shuff,'k*')
xlabel('cc measured')
ylabel('cc measured shuffled')

figure
plot(trials,measured1,'b')
hold on
plot(trials,measured2,'r')
plot(trials,signal1,'b')
plot(trials,signal2,'r')

figure
plot(res1,res2,'k*')


