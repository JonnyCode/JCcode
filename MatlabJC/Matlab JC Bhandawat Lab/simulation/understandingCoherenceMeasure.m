% understanding coherence measurments in an LN model

% JC 5/8/12

% parameters
time = [0:.0001:10] ; % sec

stim_mean = 0 ;
stim_std = 10 ;

tpeak = .05 ;   % time of peak
peakRise = .010 ;   % rise of peak
ttrough = .11 ;   % time of trough
troughDecay = .06 ;    % decay of trough
peakAmp = 1 ;    % sign and amplitude of peak
troughAmp = .2 ;   % sign and amplitude of trough

% from parameters
numPnts = length(time) ;

% stim
stim = normrnd(stim_mean,stim_std,1,numPnts) ;

% filter
Filter = simFilter(time,tpeak,peakRise,peakAmp,ttrough,troughDecay,troughAmp) ;

% linear response
RespL = conv(stim,Filter) ;
RespL = RespL(1:numPnts) ;

% Nonlinear threshold response
RespLN = RespL ;
RespLN(RespL<0)=0 ;

% coherence 
ACpStim = fft(stim).*conj(fft(stim)) ;

CCpL = fft(stim).*conj(fft(RespL)) ;
ACpRespL = fft(RespL).*conj(fft(RespL)) ;
CohL = (CCpL.*conj(CCpL))./(ACpStim.*ACpRespL) ;

CCpLN = fft(stim).*conj(fft(RespLN)) ;
ACpRespLN = fft(RespLN).*conj(fft(RespLN)) ;
CohLN = (CCpLN.*conj(CCpLN))./(ACpStim.*ACpRespLN) ;















