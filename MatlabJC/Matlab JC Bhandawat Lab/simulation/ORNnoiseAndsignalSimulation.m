% simulation to understand signal and noise separation in ORN-PN transfrom

% parameters
clear all; close all
StimType = 'pulse' ; % odor stim type (pulse, noise, filtered noise)

StimPulseTime = 0.0001 ; %(sec) zero=white noise, else square pulse of given duration

StimMean = 0 ;
StimStd = .5 ; 
NoiseMean = 0 ;
NoiseStd = .15 ;

PulseStart = 3 ; % sec
PulseStop = 3.5 ;

StimfTpeak = .150 ; % time of peak (sec)
StimfPeakRise = .050 ;  % rise of peak (std of sec)
StimfPeakAmp = 1 ; % stim values
StimfTrough = .2; % time of trough (sec)
StimfTroughDecay = .2; % decay of trough (std of sec)
StimfTroughAmp = 0; % stim values

SfTpeak = .05 ; % time of peak (sec)
SfPeakRise = .02 ;  % rise of peak (std of sec)
SfPeakAmp = 1 ; % stim values
SfTrough = .2; % time of trough (sec)
SfTroughDecay = .2; % decay of trough (std of sec)
SfTroughAmp = .05; % stim values

NfTpeak = .0001 ; % time of peak (sec)
NfPeakRise = .0001 ;  % rise of peak (std of sec)
NfPeakAmp = 1 ; % stim values
NfTtrough = .0002; % time of trough (sec)
NfTroughDecay = .0001; % decay of trough (std of sec)
NfTroughAmp = 1; % stim values

AbsRefTime = .0002 ;
RefTau = .0001 ; % refractory period decay(sec)
RefAmp = 10 ;
SpikeThresh = .5 ; % relative to standard deviations of stim and noise parameters above

SampleRate = 10000 ; % hz
StimTime = 10 ; % sec

% from params
Time = [1/SampleRate:1/SampleRate:StimTime] ;

RefPeriod = -RefAmp*exp(-Time/RefTau) ;
RefPeriod(1:AbsRefTime*SampleRate) = -inf ;

StaPnts = round((SfTrough+2*SfTroughDecay)*SampleRate) ;

% stimulus
StimFilter = simFilter(Time,StimfTpeak,StimfPeakRise,StimfPeakAmp,StimfTrough,StimfTroughDecay,StimfTroughAmp) ;

if strcmp(StimType,'noise') ;
    Stim = normrnd(0,StimStd,1,StimTime*SampleRate) ;
elseif strcmp(StimType,'filtered noise') ;
    Stim = normrnd(0,StimStd,1,StimTime*SampleRate) ;
    Stim = conv(Stim,StimFilter) ;
    Stim = Stim(1:length(Time)) ;
elseif strcmp(StimType,'pulse') ;
    Temp = zeros(1,StaPnts+StimPulseTime*SampleRate) ;
    Temp(StaPnts:StaPnts+StimPulseTime*SampleRate) = StimStd ;
    Stim = repmat(Temp,1,ceil((StimTime*SampleRate)/length(Temp))) ;
    Stim = Stim(1:length(Time)) ;
elseif strcmp(StimType,'filtered pulse') ;
    Temp = zeros(1,StaPnts+StimPulseTime*SampleRate) ;
    Temp(StaPnts:StaPnts+StimPulseTime*SampleRate) = StimStd ;
    Stim = repmat(Temp,1,ceil((StimTime*SampleRate)/length(Temp))) ;
    Stim = Stim(1:length(Time)) ;
    Stim = conv(Stim,StimFilter) ;
    Stim = Stim(1:length(Time)) ;
else
end

% noise
Noise = normrnd(0,NoiseStd,1,StimTime*SampleRate) ;

% ORN filters
SignalFilter = simFilter(Time,SfTpeak,SfPeakRise,SfPeakAmp,SfTrough,SfTroughDecay,SfTroughAmp) ;
NoiseFilter = simFilter(Time,NfTpeak,NfPeakRise,NfPeakAmp,NfTtrough,NfTroughDecay,NfTroughAmp) ;

% ORN tranduction response (generator signal)
OrnGenSignal = conv(Stim,SignalFilter) ;
OrnGenNoise = conv(Noise,NoiseFilter) ;

OrnGenSignal = OrnGenSignal(1:length(Time)) ;
OrnGenNoise = OrnGenNoise(1:length(Time)) ;

if strcmp(StimType,'noise') || strcmp(StimType,'filtered noise') ;
    OrnGenSignal = OrnGenSignal*StimStd/std(OrnGenSignal) ;
elseif strcmp(StimType,'pulse') || strcmp(StimType,'filtered pulse') ;
    OrnGenSignal = OrnGenSignal*StimStd/max(OrnGenSignal) ;
end
OrnGenNoise = OrnGenNoise*NoiseStd/std(OrnGenNoise) ; 

OrnGen = OrnGenSignal + OrnGenNoise ;

% ORN spike response (single threshold, refractory exponential)
OrnStSignal = zeros(1,length(Time)) ;
OrnStNoise = zeros(1,length(Time)) ;
OrnSt = zeros(1,length(Time)) ;

OrnGen2Signal = OrnGenSignal ;
OrnGen2Noise = OrnGenNoise ; 
OrnGen2 = OrnGen ;

for a = 1:length(Time) ;
    
    if OrnGen2Signal(a)>SpikeThresh 
        OrnStSignal(a) = 1 ;
        OrnGen2Signal(a:end) = OrnGen2Signal(a:end) + RefPeriod(1:length(Time)-a+1) ;
    end
    
    if OrnGen2Noise(a)>SpikeThresh ;
        OrnStNoise(a) = 1 ;
        OrnGen2Noise(a:end) = OrnGen2Noise(a:end) + RefPeriod(1:length(Time)-a+1) ;
    end
    
    if OrnGen2(a)>SpikeThresh ;
        OrnSt(a) = 1 ;
        OrnGen2(a:end) = OrnGen2(a:end) + RefPeriod(1:length(Time)-a+1) ;
    end
end

% power spectra

[OrnGenSignalPowerX,OrnGenSignalPower] = PowerSpectrumFinder(OrnGenSignal,SampleRate) ;
[OrnGenNoisePowerX,OrnGenNoisePower] = PowerSpectrumFinder(OrnGenNoise,SampleRate) ;
[OrnGenPowerX,OrnGenPower] = PowerSpectrumFinder(OrnGen,SampleRate) ;

[OrnStSignalPowerX,OrnStSignalPower] = PowerSpectrumFinder(OrnStSignal,SampleRate) ;
[OrnStNoisePowerX,OrnStNoisePower] = PowerSpectrumFinder(OrnStNoise,SampleRate) ;
[OrnStPowerX,OrnStPower] = PowerSpectrumFinder(OrnSt,SampleRate) ;

% estimate ORN filter without noise using STA
OrnStaSignal = staFinder(Stim,OrnStSignal,StaPnts) ;
OrnSta = staFinder(Stim,OrnSt,StaPnts) ; 

% figure

figure
plot(Time,StimFilter,'k')
hold on
plot(Time,SignalFilter,'b') ;
plot(Time,NoiseFilter,'r') ;
plot(Time(1:StaPnts),fliplr(OrnStaSignal/max(OrnStaSignal)),'b:')
plot(Time(1:StaPnts),fliplr(OrnSta/max(OrnSta)),'r:')

figure
plot(Time,OrnGenSignal,'b')
hold on
plot(Time,OrnGenNoise,'r')
plot(Time,OrnGen,'g')
plot(Time,RefPeriod,'k')
plot(Time,ones(1,length(Time))*SpikeThresh,'k:')

figure
plot(Time,OrnStSignal) 
hold on
plot(Time,OrnStNoise,'r')
plot(Time,OrnSt,'g')

figure
loglog(OrnGenSignalPowerX,OrnGenSignalPower) 
hold on
loglog(OrnGenNoisePowerX,OrnGenNoisePower,'r') 
loglog(OrnGenPowerX,OrnGenPower,'g') 

figure
loglog(OrnStSignalPowerX,OrnStSignalPower) 
hold on
loglog(OrnStNoisePowerX,OrnStNoisePower,'r') 
loglog(OrnStPowerX,OrnStPower,'g') 




