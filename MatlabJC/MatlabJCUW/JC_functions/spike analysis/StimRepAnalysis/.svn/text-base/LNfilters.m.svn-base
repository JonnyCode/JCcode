function ForIgor = LNfilters(Input,Parameters,id,A) ;

% this function will implement a linear-nonlinear model from light stimuli
% (or implied light stimuli from dynamic clamp) to spike output.

% JC 11/2/08

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

set(f1,'position',[560 420 560 420])
set(f2,'position',[0 0 560 420])
set(f3,'position',[560 0 560 420])
set(f4,'position',[0 420 560 420])
set(f5,'position',[1120 420 560 420])

% get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

% get spike data and sampling interval
epochs = str2num(Input(A).(id)) ;
for a = 1:length(epochs) ; % for each epoch  
    [predata{a}, error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get cell attached data
    [prestm{a}, error] = ITCReadEpochStm(epochs(a), 0, fp);  % get the light stimulus (this will be ignored if not cell attached)
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp);
    SI(a) = SI(a) * 1e-6; % Sampling interval in sec
end

if ~strcmp(id,'CA') % if its cell attached
% get light stimuli used to record g
cd ~/Data_analysis/Index/
if Parameters.DClight == 1 ; % if it is the first set of conductances
    cd ~/Data_analysis/Index/
    load DCGstm1 % stimuli used when recorded conductances used in dc
    prestm = DCGstm1 ; % stm
    prestm = repmat(prestm,1,ceil(length(predata)/5)) ; %repeat in a matrix as many times as needed
end
end

% interpolate all signals so they are pseudo sampled at the same rate
MSI = min(SI) ; % find the signal with the highest sampling rate
SLsec = ((Parameters.PrePnts+Parameters.StmPnts+Parameters.PostPnts)/10000) ; % entire stimulus length in sec
time = [MSI:MSI:SLsec] ; % time vector in sec 
for a = 1:length(epochs) ; % for each epoch
    data(a,:) = interp1([SI(a):SI(a):SLsec],predata{a},time,'linear','extrap') ;  % interpolate the data
    stm(a,:) = interp1([.0001:.0001:SLsec],prestm{a},time,'linear','extrap') ;  % interpolate the stimulus assuming an original sampling interval of .0001sec
end

% change points to be appropriate for new sample rate
PrePnts = round(Parameters.PrePnts*.0001/MSI) ;
StmPnts = round(Parameters.StmPnts*.0001/MSI) ;
PostPnts = round(Parameters.PostPnts*.0001/MSI) ;
STAPnts = round(Parameters.STAPnts*.0001/MSI) ;
DecimatePnts = round(Parameters.DecimatePnts*.0001/MSI) ;
ChopPnts = round(Parameters.ChopPnts*.0001/MSI) ;

% rectifier and zero stim
negstmPnts = find(stm<0) ;  % find indicies which would be getting a negative stim voltage
stm(negstmPnts) = 0 ;       % make these points zero
stm = stm - repmat(mean(stm(:,PrePnts+STAPnts:StmPnts),2),1,size(stm,2)) ; % subtract off mean of stimulus during time varying stimulus

% detect spikes in whole cell attached data
if ~strcmp(id,'CA') % if its cell attached
    SpikePnts = SpikeDetection(data,15,(1/MSI)) ; % data,threshold,samplerate
else
    SpikePnts = SpikeDetection_WC(data,-20,(1/MSI)) ; % data,threshold,samplerate
end

% get spiketrain
[PSTH,SpikeTrain] = PSTHmaker(SpikePnts,Parameters.SmoothPnts,length(data),1/MSI) ;   % PSTHmaker(SpikePnts,SmoothPnts,RawDataLength,SampleRate)

% get prespike stimuli (STA)
epochSTA = ifft(fft(SpikeTrain(:,PrePnts+STAPnts:PrePnts+StmPnts),[],2).*conj(fft(stm(:,PrePnts+STAPnts:PrePnts+StmPnts),[],2)),[],2) ; % sta is the correlation of spike train and stim
STA = sum(epochSTA,1)/sum(sum(SpikeTrain(:,PrePnts+STAPnts:PrePnts+StmPnts))) ; % sum over all the epochs and divide by the total number of spikes

figure(f1)
plot(-[MSI:MSI:length(STA)*MSI],STA)
title('STA')

%prep stim and response for lin filter
cutfactor = rem((StmPnts-STAPnts+1),ChopPnts) ; % how many pnts remain when we attempt to evenly divide the signal 
reshfactor = floor((StmPnts-STAPnts+1)/ChopPnts)*size(stm,1) ; % number of rows will we have when we divide up the signal further 

signal = stm(:,PrePnts+STAPnts:PrePnts+StmPnts-cutfactor) ; %take only the time varying stimuli and leave off any possible remainder before reshaping rows 
response = SpikeTrain(:,PrePnts+STAPnts:PrePnts+StmPnts-cutfactor) ;

signal = reshape(signal',ChopPnts,reshfactor)' ; % cut stimuli in rows into more rows
response = reshape(response',ChopPnts,reshfactor)' ;

% signal = stm(:,PrePnts+STAPnts:PrePnts+StmPnts) ; if not choping up
% response = SpikeTrain(:,PrePnts+STAPnts:PrePnts+StmPnts) ;

% get linear filter from light to spiketrain
LinearFilter = LinFilterFinder(signal,response, 1/MSI, 60) ; %(signal,response, samplerate, freqcutoff)

figure(f2)
plot(-[MSI:MSI:length(LinearFilter)*MSI],LinearFilter)
title('LinearFilter')

% get transfer function 
[powerspec_xvalues, TF] = PowerSpectrumFinder(LinearFilter,1/MSI) ; % (signal,samplerate)
TFN = TF/sum(TF) ; %  normalize by total variance
TFNS = SmoothPowerSpectrum(TFN, powerspec_xvalues, 2, 1) ; % smooth powerspec for log scaling(mean_powerspec, OriginalFreq, SmoothFact, SkipPts)

figure(f5)
plot(powerspec_xvalues, TFN)
hold on
plot(TFNS.Freq,TFNS.PowerSpec,'r')
title(['Transferfunction',id,num2str(A)])
a = gca ;
set(a,'xScale','log')
set(a,'yScale','log')

% get generator signal by convolving linear filter w/ stimulus 
LF = zeros(1,StmPnts-STAPnts) ; % prep matrix for convolution
LF(1:length(LinearFilter)) = LinearFilter ;

generator = ifft(fft(stm(:,PrePnts+STAPnts+1:PrePnts+StmPnts),[],2).*repmat(fft(LF),size(stm,1),1),[],2) ; % convolve filter with stimulus

generator = generator(:,1:StmPnts-STAPnts) ; % the end of the generator signal is an ambiguity because we are now getting the correlation with less and less data (barry's advice is to chop off end)

generator = real(generator) ;

% average spike rate for each particular range of generator values
generator_bins = min(min(generator)):.05*(max(max(generator))-min(min(generator))): max(max(generator)) ; % this makes 10 bins (11 bin points) from the values in the entire generator signal matrix
SpikeTrain2 = SpikeTrain(:,PrePnts+STAPnts:PrePnts+StmPnts) ; % take only the stimulus points

for a = 1:20 ;                % for each unique value of the generator signal
    binPnts = find(generator>generator_bins(a) & generator < generator_bins(a+1)) ;      % find all the indicies which have values within this bin
    SpksPerPnts(a) = sum(SpikeTrain2(binPnts))/length(binPnts) ;                     % find the number of spikes per sample point for each generator bin
end

figure(f3)
plot((generator(2,:)-mean(generator(2,:)))/max(generator(2,:)-mean(generator(2,:))))
hold on
plot(mean(PSTH(2:10:end,PrePnts+STAPnts:PrePnts+StmPnts))/max(mean(PSTH(2:10:end,PrePnts+STAPnts:PrePnts+StmPnts))),'r')
hold off
legend('generator','PSTH')

figure(f4)
plot(generator_bins(1:20),SpksPerPnts)
title('nonlinearityDCnoInh')

% create vectors 
identifier = ['STA',id,num2str(A)] ;
ForIgor.(identifier) = STA ;

identifier = ['LinFilter',id,num2str(A)] ;
ForIgor.(identifier) = LinearFilter ;

identifier = ['generatorBins',id,num2str(A)] ;
ForIgor.(identifier) = generator_bins(1:20) ;

identifier = ['SpksPerPnts',id,num2str(A)] ;
ForIgor.(identifier) = SpksPerPnts ;

identifier = ['MSI',id,num2str(A)] ;
ForIgor.(identifier) = MSI ;

identifier = ['TFunN',id,num2str(A)] ;
ForIgor.(identifier) = TFN ;

identifier = ['xTFun',id,num2str(A)] ;
ForIgor.(identifier) = powerspec_xvalues ;

identifier = ['TFunNS',id,num2str(A)] ;
ForIgor.(identifier) = TFNS.PowerSpec ;

identifier = ['xTFunNS',id,num2str(A)] ;
ForIgor.(identifier) = TFNS.Freq ;

end
