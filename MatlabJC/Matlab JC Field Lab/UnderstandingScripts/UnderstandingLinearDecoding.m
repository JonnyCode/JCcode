% This simulation will test and explore optimal linear decoding

% JC 5/24/18

% STIMULUS

% stim and data collection parameters (parameters of experiment) 
sampleRate = 100 ; % hz
MaxTime = 3000 ; % sec

% time 
time = [1/sampleRate:1/sampleRate:MaxTime] ;

% stimulus - This is the stimulus you are probing the system with, I provided gaussian noise but you can replace this
% function with whatever you want to probe the system with and see how that helps or hurts your ability to detect the Filter

Stim = normrnd(0,1,1,length(time)) ; % gaussian noise

% MODEL

% model parameters (parameters of biology)
tpeak = 0.05 ; % sec (gaussian filtered 0 is centered, acausal, +X becomes more causal)
peakRise = .02 ; % sec

tpeak2 = 0.02 ; % sec (gaussian filtered 0 is centered, acausal, +X becomes more causal)
peakRise2 = .01 ; % sec

% Filter - This is the filter the we are pretending the biology implements (you can replace the gaussian with any function you want). 
FilterTime = time-(MaxTime/2) ;
Filter = exp(-((FilterTime-tpeak).^2)/(2*peakRise^2)) ; % gaussian

Filter2 = -exp(-((FilterTime-tpeak2).^2)/(2*peakRise2^2)) ; % inverted gaussian

% linear response - We are pretending the system is perfectly linear at this point, so the response is just the convolution of the filter and stimulus
Response = conv(Stim,Filter,'same') ;
Spikes = Response>poissrnd(Response) ;

Response2 = conv(Stim,Filter2,'same') ;
Spikes2 = Response2>poissrnd(Response2) ;

% ANALYSIS

% bin spikes and stim
binSize = 11; % number of bins
Temp = smooth(Spikes,binSize) ;
SpikesBinned = Temp(ceil(binSize/2):binSize:end) ;
Temp = smooth(Spikes2,binSize) ;
SpikesBinned2 = Temp(ceil(binSize/2):binSize:end) ;

Temp = smooth(Stim,binSize) ;
StimBinned = Temp(ceil(binSize/2):binSize:end) ;

% try to predict stim from response using 
StimVector = Stim ;
ResponseMat = [Spikes',Spikes2'] ;
FilterLength = 50 ;

[LinFilters,Constant] = MultiCellLinFilterFinder(StimVector,ResponseMat,FilterLength) ;

StimEstimate = MultiCellLinFilterTester(LinFilters,Constant,ResponseMat) ;

StimFiltered = lowPassFilter(Stim,sampleRate,20) ;

% figure
figure
plot(Stim)
hold on
plot(StimEstimate)
%plot(StimFiltered)


