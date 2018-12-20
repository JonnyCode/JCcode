sampleNumber = 10000 ;

meanReal_range = [0:10] ;
varReal_range = [.1:1:10.1] ;

cc_range = [0:.1:1] ;

% 2-d case with correction
numTrials = length(meanReal_range)^2*length(varReal_range)*length(cc_range) ;


mean1Real = nan(1,numTrials) ;
var1Real = nan(1,numTrials) ;
mean2Real = nan(1,numTrials) ;
var2Real = nan(1,numTrials) ; 
covReal = nan(1,numTrials) ;
mean1Measured = nan(1,numTrials) ;
var1Measured = nan(1,numTrials) ;
mean2Measured = nan(1,numTrials) ;
var2Measured = nan(1,numTrials) ; 
covMeasured = nan(1,numTrials) ;
mean1Corrected = nan(1,numTrials) ;
var1Corrected = nan(1,numTrials) ;
mean2Corrected = nan(1,numTrials) ;
var2Corrected = nan(1,numTrials) ; 
covCorrected = nan(1,numTrials) ;

trial = 0 ;
for a=1:length(meanReal_range) ;
    for b=1:length(meanReal_range) ;
        for c=1:length(varReal_range) ;
            for d=1:length(cc_range) ;
                
                covSet = cc_range(d)*varReal_range(c) ;
                
                trial = trial + 1 ;
            
                % full distributions
                dist1Real = normrnd(meanReal_range(a),sqrt(varReal_range(c)),1,sampleNumber) ;
                
                dist2Real_pre = normrnd(meanReal_range(b),sqrt(varReal_range(c)),1,sampleNumber) ;
                dist2Real = dist2Real_pre*(1-(covSet/varReal_range(c)))+dist1Real*(covSet/varReal_range(c)) ;

                % rectified distributions
                dist1Measured = dist1Real ;
                dist1Measured(dist1Real<0)=0 ;

                dist2Measured = dist2Real ;
                dist2Measured(dist2Real<0)=0 ;

                % variance and covariance of full distributions
                temp = cov(dist1Real,dist2Real) ;

                mean1Real(trial) = mean(dist1Real) ;
                var1Real(trial) = temp(1,1) ;

                mean2Real(trial) = mean(dist2Real) ;
                var2Real(trial) = temp(2,2) ; 

                covReal(trial) = temp(1,2) ;
                
                % variance and covariance of rectified distributions
                temp = cov(dist1Measured,dist2Measured) ;

                mean1Measured(trial) = mean(dist1Measured) ;
                var1Measured(trial) = temp(1,1) ;

                mean2Measured(trial) = mean(dist2Measured) ;
                var2Measured(trial) = temp(2,2) ; 

                covMeasured(trial) = temp(1,2) ;


                % corrected mean and variance of rectified distributions
                [CovarRectDist, MeansRectDist] = CorrectCovarForRect(temp, [mean1Measured(trial),mean2Measured(trial)], [0,0]) ;
                
                mean1Corrected(trial) = MeansRectDist(1) ;
                var1Corrected(trial) = CovarRectDist(1,1) ;

                mean2Corrected(trial) = MeansRectDist(2) ;
                var2Corrected(trial) = CovarRectDist(2,2) ; 

                covCorrected(trial) = CovarRectDist(1,2) ;
            end
        end           
    end
end

figure
plot(mean1Real,mean1Real,'k')
hold on
plot(mean1Real,mean1Measured,'b*')
plot(mean1Real,mean1Corrected,'r*')

figure
plot(var1Real,var1Real,'k')
hold on
plot(var1Real,var1Measured,'b*')
plot(var1Real,var1Corrected,'r*')

figure
plot(covReal,covReal,'k')
hold on
plot(covReal,covMeasured,'b*')
plot(covReal,covCorrected,'r*')

