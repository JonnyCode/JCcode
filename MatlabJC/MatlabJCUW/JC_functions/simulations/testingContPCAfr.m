
% spiking pca

filterLength = 300 ;
stimulusLength = 100000;
numTrials = 100;
Thresh = 0;

spiketriggerdcov = zeros(filterLength+1,filterLength+1) ;
stimcov = zeros(filterLength+1,filterLength+1) ;

Stmfilter = simFilter([1:filterLength],10,5,10,30,10,0) ;  % simFilter(time,tpeak,peakRise,peakAmp,ttrough,troughDecay,troughAmp)

filter1 = simFilter([1:filterLength],50,30,10,100,50,1) ;  % simFilter(time,tpeak,peakRise,peakAmp,ttrough,troughDecay,troughAmp)
filter1 = filter1/sum(filter1) ;

% grahm-schmit orthogonilization
filter1diff = [diff(filter1),0] ; % derivative of filter 1 (just some way of getting another filter with proper kinetics)
filter2 = filter1diff - filter1*((filter1*filter1diff')/(filter1*filter1')) ;
filter2 = filter2/sum(filter2) ;

filter2diff = [diff(filter2),0] ;
filter3 = filter2diff - filter1*((filter1*filter2diff')/(filter1*filter1')) - filter2*((filter2*filter2diff')/(filter2*filter2')) ;
filter3 = filter3/sum(filter3) ;

figure(1);
plot([1:length(filter1)], filter1, [1:length(filter2)], filter2, [1:length(filter3)], filter3);

for trial = 1:numTrials ;

    % linear filter and stimulus
    stimulus = normrnd(0,10,1,stimulusLength) ;
    %stimulus = conv(stimulus, Stmfilter);
    stimulus = stimulus(1:stimulusLength);
    
    % generator signal
    generator1 = conv(stimulus,filter1) ;
    generator2 = conv(stimulus,filter2) ; 
    generator3 = conv(stimulus,filter3) ;

    generator1(find(generator1 < Thresh)) = 0;
    generator2(find(generator2 < Thresh)) = 0;
    generator3(find(generator3 < Thresh)) = 0;
    generator = generator1 + generator2 + generator3 ; % response is linear combination of filters

    generator = generator(1:length(stimulus)) ;
    generator(generator<0) = 0 ;
    generator = generator/max(generator) ; % non linearity

    % gererate spikes
    spikes = poissrnd(generator) ;
    spikes(spikes>1) = 1 ; 

%     figure
%     plot(generator)
%     hold on
%     plot(spikes,'r')

    % get spike triggered stim and non spike triggered stim
    spikes(1:filterLength) = 0 ; % get rid of spikes you could not use to calculate filter
    spikePnts = find(spikes==1) ; 

    spiketriggeredstimuli = nans(length(spikePnts),filterLength+1) ;
    arbitrarystimuli = nans(length(spikePnts),filterLength+1) ;
    spiketriggeredRes = nans(length(spikePnts),filterLength+1) ;

    for a=1:length(spikePnts) ; % for every spike 
       spiketriggeredstimuli(a,:) = stimulus(spikePnts(a)-filterLength:spikePnts(a)) ;
       randomPnt = randi(length(spikePnts)-filterLength)+filterLength ; 
       arbitrarystimuli(a,:) = stimulus(randomPnt-filterLength:randomPnt) ;
    end

    % sta
    sta = mean(spiketriggeredstimuli) ;
    arbitraryAverage = mean(arbitrarystimuli) ;

%     figure
%     plot(sta)
%     hold on
%     plot(fliplr(filter1),'r')
%     plot(fliplr(filter2),'y')
%     plot(fliplr(filter3),'k')
%     plot(arbitraryAverage,'g')


    % covariance
    for a=1:length(spikePnts) ; % for every spike 
        spiketriggeredRes(a,:) = spiketriggeredstimuli(a,:)-sta ;
    end

    spiketriggerdcov = spiketriggerdcov + cov(spiketriggeredRes) ;
    stimcov = stimcov + cov(arbitrarystimuli) ;

end
  
spiketriggerdcov = spiketriggerdcov/numTrials ;
stimcov = stimcov/numTrials ;

% figure;
% mesh(spiketriggerdcov);
% 
% figure;
% mesh(stimcov);


figure(2) ;
mesh(spiketriggerdcov-stimcov) ;

% eigenmodes
[temp, EigVal] = eig(spiketriggerdcov-stimcov);
[EigVec, temp] = eigs(spiketriggerdcov-stimcov, 5);

Indices = find(EigVal(:) > 0);

figure(3)
%plot(sort(abs(EigVal(Indices)), 'descend'), 'o'); % freds
plot(sort(abs(EigVal), 'descend'), 'o'); % mine

figure(4); clf
plot(EigVec);

%%



% continuous PCA

filterLength = 300 ;
stimulusLength = 100000;
numTrials = 100;
Thresh = 0;

spiketriggerdcov = zeros(filterLength+1,filterLength+1) ;
stimcov = zeros(filterLength+1,filterLength+1) ;

Stmfilter = simFilter([1:filterLength],10,5,10,30,10,0) ;  % simFilter(time,tpeak,peakRise,peakAmp,ttrough,troughDecay,troughAmp)

filter1 = simFilter([1:filterLength],50,30,10,100,50,1) ;  % simFilter(time,tpeak,peakRise,peakAmp,ttrough,troughDecay,troughAmp)
filter1 = filter1/sum(filter1) ;

% grahm-schmit orthogonilization
filter1diff = [diff(filter1),0] ; % derivative of filter 1 (just some way of getting another filter with proper kinetics)
filter2 = filter1diff - filter1*((filter1*filter1diff')/(filter1*filter1')) ;
filter2 = filter2/sum(filter2) ;

filter2diff = [diff(filter2),0] ;
filter3 = filter2diff - filter1*((filter1*filter2diff')/(filter1*filter1')) - filter2*((filter2*filter2diff')/(filter2*filter2')) ;
filter3 = filter3/sum(filter3) ;

figure(1);
plot([1:length(filter1)], filter1, [1:length(filter2)], filter2, [1:length(filter3)], filter3);

Final_ResponseTrigCov = zeros(filterLength,filterLength) ;

for trial = 1:numTrials ;

    % linear filter and stimulus
    stimulus = normrnd(0,10,1,stimulusLength) ;
    %stimulus = conv(stimulus, Stmfilter);
    stimulus = stimulus(1:stimulusLength);
    
    % generator signal
    generator1 = conv(stimulus,filter1) ;
    generator2 = conv(stimulus,filter2) ; 
    generator3 = conv(stimulus,filter3) ;

    generator1(find(generator1 < Thresh)) = 0;
    generator2(find(generator2 < Thresh)) = 0;
    generator3(find(generator3 < Thresh)) = 0;
    generator = generator1 + generator2 + generator3 ; % response is linear combination of filters

    generator = generator(1:length(stimulus)) ;
    generator(generator<0) = 0 ;
    generator = generator/max(generator) ; % non linearity

    % covariance
    genMean = mean(generator) ;
    
    STA = zeros(1,filterLength) ;
    ResponseTrigCov = zeros(filterLength,filterLength) ;
    StimCov = zeros(filterLength,filterLength) ;
    
    for t=filterLength/2+1:stimulusLength-filterLength/2 ;
        for Tau = -filterLength/2:filterLength/2-1 ;
            
            STA(Tau+filterLength/2+1) = STA(Tau+filterLength/2+1) + generator(t)*stimulus(t-Tau) ;
            
            for Tau2 = -filterLength/2:filterLength/2-1 ;   
                ResponseTrigCov(Tau+filterLength/2+1,Tau2+filterLength/2+1) = ResponseTrigCov(Tau+filterLength/2+1,Tau2+filterLength/2+1) + (generator(t))*stimulus(t-Tau)*stimulus(t-Tau2) ;
                StimCov(Tau+filterLength/2+1,Tau2+filterLength/2+1) = StimCov(Tau+filterLength/2+1,Tau2+filterLength/2+1) + (genMean)*stimulus(t-Tau)*stimulus(t-Tau2) ;
            end
        end
    end

    STA = STA/length([filterLength/2+1:stimulusLength-filterLength/2]) ;
    ResponseTrigCov = ResponseTrigCov/length([filterLength/2+1:stimulusLength-filterLength/2]) ;
    StimCov = StimCov/length([filterLength/2+1:stimulusLength-filterLength/2]) ;

    for Tau = -filterLength/2:filterLength/2-1 ; 
        for Tau2 = -filterLength/2:filterLength/2-1 ; 
            ResponseTrigCov(Tau+filterLength/2+1,Tau2+filterLength/2+1) = ResponseTrigCov(Tau+filterLength/2+1,Tau2+filterLength/2+1)....
                - STA(Tau+filterLength/2+1)*STA(Tau2+filterLength/2+1) - StimCov(Tau+filterLength/2+1,Tau2+filterLength/2+1) ;
        end
    end
    
            
    Final_ResponseTrigCov = Final_ResponseTrigCov + ResponseTrigCov ;

end
  
Final_ResponseTrigCov = Final_ResponseTrigCov/numTrials ;

% eigenmodes
[temp, EigVal] = eig(Final_ResponseTrigCov);
[EigVec, temp] = eigs(Final_ResponseTrigCov, 5);

Indices = find(EigVal(:) > 0);

figure(3)
%plot(sort(abs(EigVal(Indices)), 'descend'), 'o'); % freds
plot(sort(abs(EigVal), 'descend'), 'o'); % mine

figure(4); clf
plot(EigVec);









