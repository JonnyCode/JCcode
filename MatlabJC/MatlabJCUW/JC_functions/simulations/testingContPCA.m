
% spiking pca

filterLength = 200 ;
stimulusLength = 1000000 ;
numTrials = 100 ;

spiketriggerdcov = zeros(filterLength+1,filterLength+1) ;
stimcov = zeros(filterLength+1,filterLength+1) ;

for trial = 1:numTrials ;

    % linear filter and stimulus
    stimulus = normrnd(0,10,1,stimulusLength) ;

    filter1 = simFilter([1:filterLength],50,30,10,100,50,1) ;  % simFilter(time,tpeak,peakRise,peakAmp,ttrough,troughDecay,troughAmp)
    filter1 = filter1/sum(filter1) ;

    % grahm-schmit orthogonilization
    filter1diff = [diff(filter1),0] ; % derivative of filter 1 (just some way of getting another filter with proper kinetics)
    filter2 = filter1diff - filter1*((filter1*filter1diff')/(filter1*filter1')) ;
    filter2 = filter2/sum(filter2) ;

    filter2diff = [diff(filter2),0] ;
    filter3 = filter2diff - filter1*((filter1*filter2diff')/(filter1*filter1')) - filter2*((filter2*filter2diff')/(filter2*filter2')) ;
    filter3 = filter3/sum(filter3) ;

    % generator signal
    generator1 = conv(stimulus,filter1) ;
    generator2 = conv(stimulus,filter2) ; 
    generator3 = conv(stimulus,filter3) ;
    %generator2 = 0 ;
    %generator3 = 0 ;
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


figure ;
mesh(spiketriggerdcov-stimcov) ;

% eigenmodes
[EigVec, EigVal] = eig(spiketriggerdcov-stimcov);

figure
semilogy(abs(EigVal), 'o');

figure
plot(EigVec(:,1),'b');
hold on
plot(EigVec(:,2),'r');
plot(EigVec(:,3),'g');





% % continuous PCA
% 
% stimulus = normrnd(0,10,1,10000) ;
% 
% filter1 = simFilter([1:200],50,30,10,100,50,1) ;  % simFilter(time,tpeak,peakRise,peakAmp,ttrough,troughDecay,troughAmp)
% 
% filterNullSet = null(filter1) ; % null basis of filter 1
% filter2 = filterNullSet(:,100) ;  % arbitray chooce of a null filter from null basis
% 
% response1 = conv(stimulus,filter1) ;
% response2 = conv(stimulus,filter2) ; 
% response = response1 + response2 ; % response is linear combination of filters
% response = response(1:length(stimulus)) ;
% 
% 
% [covar,EigVal,EigVec] = contPCA(response,stimulus,200) ; %contPCA(response,stimulus,lag)





