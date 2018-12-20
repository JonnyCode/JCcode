function ForIgor = LNstmRep(Input,Parameters,id,A) ;

% this function will create use a linear non-linear model to predict the
% light stimuli given a spike output and assess signal to noise ratio 

%JC 10/15/08

% prep figure positions
figure
f1=gcf ;

figure
f2=gcf ;

figure
f3=gcf ;

figure
f4=gcf ;

figure
f5=gcf ;

figure
f6=gcf ;

set(f1,'position',[560 420 560 420])
set(f2,'position',[0 0 560 420])
set(f3,'position',[560 0 560 420])
set(f4,'position',[0 420 560 420])
set(f5,'position',[1120 0 560 420])
set(f6,'position',[1120 420 560 420])

% get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

% get spike data and sampling interval
epochs = str2num(Input(A).(id)) ;
for a = 1:length(epochs) ; % for each epoch  
    [predata{a}, error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get cell attached data
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp);
    SI(a) = SI(a) * 1e-6; % Sampling interval in sec
end

% get light stimuli used to record g
cd ~/Data_analysis/Index/
load DCGstm % stimuli used when recorded conductances used in dc
prestm = DCGstm ; % stm

% rectifier and zero stim
negstmPnts = find(prestm<0) ;  % find indicies which would be getting a negative stim voltage
prestm(negstmPnts) = 0 ;       % make these points zero
prestm = prestm - repmat(mean(prestm(:,Parameters.PrePnts:Parameters.StmPnts),2),1,size(prestm,2)) ; % subtract off mean of stimulus during time varying stimulus
prestm = repmat(prestm,ceil(length(predata)/size(prestm,1)),1) ; %repeat in a matrix as many times as needed

% interpolate all signals so they are pseudo sampled at the same rate
MSI = min(SI) ; % find the signal with the highest sampling rate
SLsec = ((Parameters.PrePnts+Parameters.StmPnts+Parameters.PostPnts)/10000) ; % entire stimulus length in sec
time = [MSI:MSI:SLsec] ; % time vector in sec 
for a = 1:length(epochs) ; % for each epoch
    data(a,:) = interp1([SI(a):SI(a):SLsec],predata{a},time,'linear','extrap') ;  % interpolate the data
    stm(a,:) = interp1([.0001:.0001:SLsec],prestm(a,:),time,'linear','extrap') ;  % interpolate the stimulus assuming an original sampling interval of .0001sec
end

% change points to be appropriate for new sample rate
PrePnts2 = round(Parameters.PrePnts*.0001/MSI) ;
StmPnts2 = round(Parameters.StmPnts*.0001/MSI) ;
PostPnts2 = round(Parameters.PostPnts*.0001/MSI) ;
STA_pnts2 = round(Parameters.STAPnts*.0001/MSI) ;
DecimatePnts2 = round(Parameters.DecimatePnts*.0001/MSI) ;

% detect spikes in whole cell attached data
SpikePnts = SpikeDetection_WC(data,-20,(1/MSI)) ; % data,threshold,samplerate

% get spiketrain
[PSTH,SpikeTrain] = PSTHmaker(SpikePnts,Parameters.SmoothPnts,length(data),1/MSI) ;   % PSTHmaker(SpikePnts,SmoothPnts,RawDataLength,SampleRate)

% get linear filter from spikes to stimulus
LinearFilter = LinFilterFinder(SpikeTrain(:,PrePnts2:StmPnts2+PrePnts2),stm(:,PrePnts2:StmPnts2+PrePnts2), 1/MSI, 1/(MSI*2)) ; %(signal,response, samplerate, freqcutoff)

% convolve linear filter with spike train
for a=1:length(epochs) ; %for each spike train
LinPredStm(a,:) = real(ifft(fft(LinearFilter).*fft(SpikeTrain(a,PrePnts2:StmPnts2+PrePnts2)))) ;
end

stm2 = stm(:,PrePnts2:StmPnts2+PrePnts2) ; % take those points during stimulus only

figure(f1)
plot([1:length(LinearFilter)],LinearFilter)
title('linear filter from spikes to light')

% find nonlinearity and make linear-nonlinear prediction
GenBin = min(min(LinPredStm)):.05*(max(max(LinPredStm))-min(min(LinPredStm))): max(max(LinPredStm)) ; % this makes 20 bins (21 bin points) from the values in the entire generator signal matrix

for a = 1:length(GenBin)-1 ;                % for each unique value of the generator signal
    if a<length(GenBin)-1
        binPnts{a} = find(LinPredStm>=GenBin(a) & LinPredStm < GenBin(a+1)) ;     % find all the indicies which have values within this bin 
    else
        binPnts{a} = find(LinPredStm>=GenBin(a) & LinPredStm <= GenBin(a+1)) ;     % this if statment makes sure we include the maximum in our last bin
    end
    numPnts(a) = length(binPnts{a}) ;
    TrueStmMean(a) = sum(stm2(binPnts{a}))/numPnts(a) ;             % find the mean value of the true stiumulus  per sample point for each generator bin
    TrueStmMed(a) = median(stm2(binPnts{a})) ;                      % find the median value of the true stim
    TrueStmMod(a) = mode(stm2(binPnts{a})) ;                        % find the mode of the true stim
    %LNPredStm(binPnts{a}) = TrueStm(a) ;                        % make the linear nonlinear prediction
end

LNPredStm = LinPredStm ; % prep the linear nonlinear prediction matrix for speed
GenBin2 = min(min(LinPredStm)):.001*(max(max(LinPredStm))-min(min(LinPredStm))): max(max(LinPredStm)) ; % bin unique values of the linear prediction


for a=1:length(GenBin2)-1 ; % for each unique value bin of the generator signal
    if a<999
        UniPnts = find(LinPredStm>=GenBin2(a) & LinPredStm<GenBin2(a+1)) ; % find the indicies that are this in this bin
    else
        UniPnts = find(LinPredStm>=GenBin2(a) & LinPredStm<=GenBin2(a+1)) ; % this if statment makes sure to include the last point in our last bin
    end
    LNPredStm(UniPnts) = interp1(GenBin(1:20),TrueStmMean,(GenBin2(a)+GenBin2(a+1))/2,'linear','extrap') ; % make the prediction based off the nonlinearity using the center of the bin
end

for a=1:5 ; % for each unique light stim
    LNPredStmMean(a,:) = mean(LNPredStm(a:5:end,:),1) ; % get the mean predicted stim
    LNPredStmVar(a,:) = var(LNPredStm(a:5:end,:),1) ; % time dependent variance of pred stim
end
    

figure(f2) 
plot(GenBin(1:20),TrueStmMean,'*-')
xlabel('GenBin')
ylabel('stm mean')
title('nonlinearity')

for a=1:size(stm2,1)
figure(f3)
plot(stm2(a,:),'b')
hold on
plot(LinPredStm(a,:),'r')
plot(LNPredStm(a,:),'g')
legend('actual','lin pred','lin-nonlin pred')
title('predictions from spike train')
hold off
end

%stimulus distributions
minBin = min(min(stm2(1:end)),min(LNPredStm(1:end))) ;
maxBin = max(max(stm2(1:end)),max(LNPredStm(1:end))) ;
diffBin = (maxBin-minBin)*.001 ;
DisBins = minBin:diffBin:maxBin+diffBin ;
predDis = hist(LNPredStm(1:end),DisBins) ;
stmDis = hist(stm2(1:end),DisBins) ;

diffDis = stmDis-predDis ;

figure(f4)
plot(DisBins,stmDis,'*')
hold on
plot(DisBins,predDis,'g*')
hold on
plot(DisBins,diffDis,'r')
hold off
legend('stm','pred','diff') 
title('light distributions')

% find residuals and signal to noise ratio of estimation
Resid = stm2 - LNPredStm ; % actual - predicted stimulus
SNR = mean(var(stm2,[],2)./var(Resid,[],2)) ; % variance of stim / variance of residual error averaged over all trials

Resid2 = stm2 - LinPredStm ;
SNR2 = mean(var(stm2,[],2)./var(Resid2,[],2)) ;

% find gain and signal to noise ratio without confound of systematic errors

binsize = floor(length(LNPredStm)/6) ; % segment size if we divide each estimate into 6 bins
bincut = rem(length(LNPredStm),6) ; % the part we should cut off so the we get even segments
LNPredStmShort = LNPredStm(:,1:end-bincut) ; % make the cut for estimate
StmShort = stm2(:,1:end-bincut) ; % make cut for actual light stim
binnum = numel(LNPredStmShort)/binsize ; % number of total segments
PredBins = reshape(LNPredStmShort', binsize, binnum)' ; % divide the estimate into segments
StmBins = reshape(StmShort', binsize, binnum)' ; % divide the light into segments

StmBinsAmp = sqrt(2*real(fft(StmBins,[],2).*conj(fft(StmBins,[],2)))) ;  % find the amplitude of each fourier component in each segment
PredBinsAmp = sqrt(2*real(fft(PredBins,[],2).*conj(fft(PredBins,[],2)))) ;

maxfreq = (1/MSI)/2 ;
powerspec_xvalues = [0:2*maxfreq/binsize:maxfreq]; % time vector for power spectrum


StmBinsAmp = StmBinsAmp(:,1:length(StmBinsAmp)/2) ; % the fft has neg frequencies (we accounted for this above by multiplying by 2)
PredBinsAmp = PredBinsAmp(:,1:length(PredBinsAmp)/2) ;

for a=1:find(powerspec_xvalues>60,1,'first') %length(StmBinsAmp2) ; % for each fourier component
    figure(f5)
    plot(StmBinsAmp(:,a),PredBinsAmp(:,a),'k.') ; % plot a point for each segment
    bEst = (min(PredBinsAmp(:,a))+max(PredBinsAmp(:,a)))/2 ; % estmate y intercept
    Coeffs = nlinfit(StmBinsAmp(:,a),PredBinsAmp(:,a),@StraightLineFit,[.5,bEst]) ; % fit straight line to data
    range = [min(StmBinsAmp(:,a)):max(StmBinsAmp(:,a))] ; % range of x values
    hold on
    plot(range,(range*Coeffs(1))+repmat(Coeffs(2),1,length(range)),'r') %plot best fit line
    neff{a} = (PredBinsAmp(:,a)-repmat(Coeffs(2),length(PredBinsAmp(:,a)),1))/Coeffs(1) - StmBinsAmp(:,a) ; % the distribution around the line of best fit (the effective noise)
    Neff(a) = var(neff{a}) ; % the variance of the noise
    S(a) = var(StmBinsAmp(:,a)) ; % the variance of the signal
    SNRfourier(a) = S(a)/Neff(a) ; % signal to noise
    hold off
end  

% integrate SNR over frequency components to find lower bound of
% information rate
FourierStep = 2*maxfreq/binsize ; % step size
InfoRate = (sum(log2(ones(1,50)+SNRfourier(1:50)))*FourierStep)/(4*pi) ;

% find maximum entropy of system
prob = sum(sum(SpikeTrain(:,PrePnts2:StmPnts2+PrePnts2),2))/(numel(SpikeTrain(:,PrePnts2:StmPnts2+PrePnts2))*(MSI*1000)) ; %probability of spike per ms
deltaT = MSI*1000;
Entropy = -((1-prob)*log2(1-prob)+prob*log2(prob))/deltaT ;

figure(f6)
plot(powerspec_xvalues(1:length(SNRfourier)),SNRfourier)
text(50,1,['Entropy=',num2str(Entropy)])
text(50,.5,['InfoRate=',num2str(InfoRate)])

% prep matricies for igor export
identifier = ['Entropy',id,num2str(A)] ;
ForIgor.(identifier) = Entropy ;

identifier = ['InfoRate',id,num2str(A)] ;
ForIgor.(identifier) = InfoRate ;

identifier = ['SNRx',id,num2str(A)] ;
ForIgor.(identifier) = powerspec_xvalues(1:length(SNRfourier)) ;

identifier = ['SNR',id,num2str(A)] ;
ForIgor.(identifier) = SNR ;

identifier = ['DisX',id,num2str(A)] ;
ForIgor.(identifier) = DisBins ;

identifier = ['stmDis',id,num2str(A)] ;
ForIgor.(identifier) = stmDis ;

identifier = ['predDis',id,num2str(A)] ;
ForIgor.(identifier) = predDis ;

identifier = ['diffDis',id,num2str(A)] ;
ForIgor.(identifier) = diffDis ;

identifier = ['PredStm',id,num2str(A)] ;
ForIgor.(identifier) = reshape(LNPredStmMean',1,numel(LNPredStmMean)) ;

identifier = ['MSI',id,num2str(A)] ;
ForIgor.(identifier) = MSI ;

end







