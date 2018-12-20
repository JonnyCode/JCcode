% simulation to understand signal and noise separation in ORN-PN transfrom.
% Modified from "ORNnoiseAndsignalSimulution"

% questions to address
% 1) How well can Orn and PN filter be estimated under different stimulus conditions?
% 2) How well does optimal LN model seperate signal from noise?
% 3) How do optimal linear and nonlinear filters change with noise mean?

clear all; close all

% time parameters
SampleRate = 10000 ; % hz
StimTime = 100 ; % sec

% stim (in terms of spike output)
StimType = 'noise' ; % odor stim type (pulse, noise, filtered noise)

PulseStart = 3 ; % sec
PulseStop = 3.5 ;
StimPulseTime = 0.0001 ; %(sec) zero=white noise, else square pulse of given duration

StimMean = 45 ; % (Hz) average spike number
StimStd = 25 ; % (Hz) std of spike number
NoiseMean = 20 ; % (Hz) noise var = noise mean

AbsRefTime = .001 ; % (sec) also sets time bin

% linear filter parameters
StimfTpeak = .150 ; % time of peak (sec)
StimfPeakRise = .050 ;  % rise of peak (std of sec)
StimfPeakAmp = 1 ; % stim values
StimfTrough = .2 ; % time of trough (sec)
StimfTroughDecay = .2 ; % decay of trough (std of sec)
StimfTroughAmp = 0 ; % stim values

SfTpeak = .05 ; % time of peak (sec)
SfPeakRise = .02 ;  % rise of peak (std of sec)
SfPeakAmp = 1 ; % stim values
SfTrough = .2 ; % time of trough (sec)
SfTroughDecay = .2 ; % decay of trough (std of sec)
SfTroughAmp = .05 ; % stim values

NfTpeak = .0001 ; % time of peak (sec)
NfPeakRise = .0001 ;  % rise of peak (std of sec)
NfPeakAmp = 1 ; % stim values
NfTtrough = .0002 ; % time of trough (sec)
NfTroughDecay = .0001; % decay of trough (std of sec)
NfTroughAmp = 1 ; % stim values

% from params
Time = [1/SampleRate:1/SampleRate:StimTime] ;

StaPnts = round((SfTrough+2*SfTroughDecay)*SampleRate) ;

% stimulus
StimFilter = simFilter(Time,StimfTpeak,StimfPeakRise,StimfPeakAmp,StimfTrough,StimfTroughDecay,StimfTroughAmp) ;

if strcmp(StimType,'noise') ;
    Stim = normrnd(0,1,1,StimTime*SampleRate) ;
elseif strcmp(StimType,'filtered noise') ;
    Stim = normrnd(0,1,1,StimTime*SampleRate) ;
    Stim = conv(Stim,StimFilter) ;
    Stim = Stim(1:length(Time)) ;
elseif strcmp(StimType,'pulse') ;
    Temp = zeros(1,StaPnts+StimPulseTime*SampleRate) ;
    Temp(StaPnts:StaPnts+StimPulseTime*SampleRate) = 1 ;
    Stim = repmat(Temp,1,ceil((StimTime*SampleRate)/length(Temp))) ;
    Stim = Stim(1:length(Time)) ;
elseif strcmp(StimType,'filtered pulse') ;
    Temp = zeros(1,StaPnts+StimPulseTime*SampleRate) ;
    Temp(StaPnts:StaPnts+StimPulseTime*SampleRate) = 1 ;
    Stim = repmat(Temp,1,ceil((StimTime*SampleRate)/length(Temp))) ;
    Stim = Stim(1:length(Time)) ;
    Stim = conv(Stim,StimFilter) ;
    Stim = Stim(1:length(Time)) ;
else
end

% noise
Noise = normrnd(0,1,1,StimTime*SampleRate) ;

% ORN filters
SignalFilter = simFilter(Time,SfTpeak,SfPeakRise,SfPeakAmp,SfTrough,SfTroughDecay,SfTroughAmp) ;
NoiseFilter = simFilter(Time,NfTpeak,NfPeakRise,NfPeakAmp,NfTtrough,NfTroughDecay,NfTroughAmp) ;

% ORN tranduction response (generator signal)
OrnGenSignal = conv(Stim,SignalFilter) ;
OrnGenNoise = conv(Noise,NoiseFilter) ;

OrnGenSignal = OrnGenSignal(1:length(Time)) ;
OrnGenNoise = OrnGenNoise(1:length(Time)) ;

% ORN generator signal into requested firing rate (as probability per time bin)
if strcmp(StimType,'noise') || strcmp(StimType,'filtered noise') ;
    OrnGenSignal = OrnGenSignal - mean(OrnGenSignal) ;
    OrnGenSignal = OrnGenSignal*((StimStd*AbsRefTime)/std(OrnGenSignal)) ; % adjust std
    OrnGenSignal = OrnGenSignal+(StimMean*AbsRefTime) ; % adjust mean
elseif strcmp(StimType,'pulse') || strcmp(StimType,'filtered pulse') ;
    OrnGenSignal = OrnGenSignal*(StimMean*AbsRefTime)/max(OrnGenSignal) ; % adjust amplitude
end

OrnGenNoise = OrnGenNoise - mean(OrnGenNoise) ;
OrnGenNoise = OrnGenNoise*(sqrt(NoiseMean)*AbsRefTime)/std(OrnGenNoise) ; % adjust std of noise
OrnGenNoise = OrnGenNoise+(NoiseMean*AbsRefTime) ; % adjust mean of noise

OrnGen = OrnGenSignal + OrnGenNoise ;

% ORN nonlinearity (linear with thresh and sat)
OrnSpikeProbSignal = OrnGenSignal ;
OrnSpikeProbNoise = OrnGenNoise ;
OrnSpikeProb = OrnGen ;

OrnSpikeProbSignal(OrnSpikeProbSignal>1)=1;
OrnSpikeProbSignal(OrnSpikeProbSignal<0)=0;

OrnSpikeProbNoise(OrnSpikeProbNoise>1)=1 ;
OrnSpikeProbNoise(OrnSpikeProbNoise<0)=0 ;

OrnSpikeProb(OrnSpikeProb>1)=1;
OrnSpikeProb(OrnSpikeProb<0)=0;

% ORN spike response 
OrnStSignal = zeros(1,length(Time)) ;
OrnStNoise = zeros(1,length(Time)) ;
OrnSt = zeros(1,length(Time)) ;

RandVec = rand(1,StimTime*SampleRate) ; % random uniform vector (0-1)

for a = round([AbsRefTime:AbsRefTime:StimTime]*SampleRate) ; % for each time bin
    
    if OrnSpikeProbSignal(a)>RandVec(a) 
        OrnStSignal(a) = 1 ;
    end
    
    if OrnSpikeProbNoise(a)>RandVec(a) ;
        OrnStNoise(a) = 1 ;
    end
    
    if OrnSpikeProb(a)>RandVec(a) ;
        OrnSt(a) = 1 ;
    end
end

% optimal linear noise filter from orn spike probability 
OrnSpikeProbSignalAc = xcorr(OrnSpikeProbSignal,SampleRate) ;
OrnSpikeProbNoiseAc = xcorr(OrnSpikeProbNoise,SampleRate) ;

OptLF = OrnSpikeProbSignalAc./(OrnSpikeProbSignalAc+OrnSpikeProbNoiseAc) ;

% Orn spike probability histograms
OrnHistX = [min([OrnSpikeProbNoise,OrnGenSignal]):range([OrnSpikeProbNoise,OrnSpikeProbSignal])/1000:max([OrnSpikeProbNoise,OrnSpikeProbSignal])] ;

OrnSpikeProbSignalHist = hist(OrnSpikeProbSignal,OrnHistX) ;
OrnSpikeProbNoiseHist = hist(OrnSpikeProbNoise,OrnHistX) ;

% optimal nonlinear noise filter
PnGenSignal = conv(OptLF,OrnSpikeProbSignal) ;
PnGenSignal = PnGenSignal(1:length(Time)) ;

PnGenNoise = conv(OptLF,OrnSpikeProbNoise) ;
PnGenNoise = PnGenNoise(1:length(Time)) ;

PnGen = conv(OptLF,OrnSpikeProb) ;
PnGen = PnGenNoise(1:length(Time)) ;

HistX = [min([PnGenNoise,PnGenSignal]):range([PnGenNoise,PnGenSignal])/1000:max([PnGenNoise,PnGenSignal])] ;

PnGenSignalHist = hist(PnGenSignal,HistX) ;
PnGenNoiseHist = hist(PnGenNoise,HistX) ;

OptNlF = PnGenSignalHist./(PnGenSignalHist+PnGenNoiseHist) ;

% PN nonlinearity
PnSpikeProbSignal = nan(1,length(Time)) ;
PnSpikeProbNoise = nan(1,length(Time)) ;
PnSpikeProb = nan(1,length(Time)) ;

for a=1:length(Time) ;
    [Min,MinI] = min(HistX-PnGenSignal(a)) ;
    PnSpikeProbSignal(a) = PnGenSignal(a)*OptNlF(MinI) ;
    
    [Min,MinI] = min(HistX-PnGenNoise(a)) ;
    PnSpikeProbNoise(a) = PnGenNoise(a)*OptNlF(MinI) ;

    [Min,MinI] = min(HistX-PnGen(a)) ;
    PnSpikeProb(a) = PnGen(a)*OptNlF(MinI) ;    
end

% Pn spike prob hist
PnSpikeProbSignalHist = hist(PnSpikeProbSignal,HistX) ;
PnSpikeProbNoiseHist = hist(PnSpikeProbNoise,HistX) ;
PnSpikeProbHist = hist(PnSpikeProb,HistX) ;


% power spectra (with lowest frequency at LowFreq Hz)
LowFreq = 5 ;

% cut up into smaller chunks
BStep = SampleRate/LowFreq ; 

for a=1:floor(length(Time)/BStep) ;
    OrnStSignalBlocks(a,:) = OrnStSignal((a-1)*BStep+1:a*BStep) ;
    OrnStNoiseBlocks(a,:) = OrnStNoise((a-1)*BStep+1:a*BStep) ;
    OrnStBlocks(a,:) = OrnSt((a-1)*BStep+1:a*BStep) ;
end

[OrnGenSignalPowerX,OrnGenSignalPower] = PowerSpectrumFinder(OrnGenSignal,SampleRate) ;
[OrnGenNoisePowerX,OrnGenNoisePower] = PowerSpectrumFinder(OrnGenNoise,SampleRate) ;
[OrnGenPowerX,OrnGenPower] = PowerSpectrumFinder(OrnGen,SampleRate) ;

[OrnStSignalPowerX,OrnStSignalPower] = PowerSpectrumFinder(OrnStSignalBlocks,SampleRate) ;
[OrnStNoisePowerX,OrnStNoisePower] = PowerSpectrumFinder(OrnStNoiseBlocks,SampleRate) ;
[OrnStPowerX,OrnStPower] = PowerSpectrumFinder(OrnStBlocks,SampleRate) ;

% estimate ORN filter using STA
OrnStaSignal = staFinder(Stim,OrnStSignal,StaPnts) ;
OrnSta = staFinder(Stim,OrnSt,StaPnts) ; 


% figure

figure
plot(Time,StimFilter,'k')
hold on
plot(Time,SignalFilter,'b') ;
plot(Time,NoiseFilter,'r') ;
plot(Time(1:StaPnts),fliplr(OrnStaSignal/max(OrnStaSignal)),'b:')
plot(Time(1:StaPnts),fliplr(OrnSta/max(OrnSta)),'g:')

figure
subplot(3,1,1)
plot(OrnHistX,OrnSpikeProbSignalHist,'b')
hold on
plot(OrnHistX,OrnSpikeProbNoiseHist,'r')

subplot(3,1,2)
plot(HistX,PnGenSignalHist,'b')
hold on
plotyy(HistX,PnGenNoiseHist,HistX,OptNlF)
plot(HistX,OptNlF,'k')

figure
plot(HistX,PnSpikeProbSignalHist)

figure
plot(Time,OrnSpikeProb)
hold on
plot(Time,PnSpikeProb,'r')


