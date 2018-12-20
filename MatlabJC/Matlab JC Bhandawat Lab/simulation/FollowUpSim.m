% correlation and convolution simulation, This simulation will explore
% finding filters and prediciting responses in a linear system
% JC 9/18/12

% parameters
sampleRate = 10000 ; % hz
MaxTime = 10 ; % sec

tpeak = .01 ; % sec (gaussian filtered 0 is centered, acausal, +X becomes more causal)
peakRise = .02 ; % sec

PulseTime = 2 ; % sec

StimType = 'Pulse' ; % 'Pulse' or nothing if you want gaussian noise
lpfCutOff = 5000; % hz , if you want to slow your probe stim down

MaxLag = 1 ; % sec, time length of filter you are finding

% time
time = [1/sampleRate:1/sampleRate:MaxTime] ;

% MODEL
% Filter - This is the filter the we are pretending the biology implements (you can replace the gaussian with any function you want). 
Filter = exp(-((time-(tpeak+round(MaxTime/2))).^2)/(2*peakRise^2)) ; % gaussian

% stimulus - This is the stimulus you are probing the system with, I provided  either a pulse or gaussian noise but you can replace these
% functions with whatever you want to probe the system with and see how that helps or hurts your ability to detect the Filter
if strcmp(StimType,'Pulse') ; % pulse
    Stim = zeros(1,length(time)) ;
    Stim(PulseTime*sampleRate) = 1 ;
else 
    Stim = normrnd(0,1,1,length(time)) ; % gaussian noise
end

% low pass filter - this will low pass the stim you are probing the system with the lpfCutOff and see how it effects your ability to detect the filter
if lpfCutOff<sampleRate/2 ; % niquist
    % implement a butterworth 
    Wn = lpfCutOff*2/sampleRate; %normalized frequency cutoff (fraction of nyquist)
    [z, p, k] = butter(1,Wn,'low'); %high can be changed to low
    [sos,g]=zp2sos(z,p,k); % convesion?
    myfilt=dfilt.df2sos(sos,g);
    Stim = filter(myfilt,Stim')'; %implementation
end

% response - We are pretending the system is perfectly linear, so the response is just the convolution of the filter and stimulus
Response = conv(Stim,Filter) ;
Response = Response(length(time)/2:end-length(time)/2) ;

% FINDING FILTER - As the expereimentor you only get the stimulus you provided and the response you got and you are trying to find the linear filter
cc= xcorr(Stim,Response,MaxLag*sampleRate) ;
ccTime = ([1:length(cc)] - (length(cc)+1)/2)/sampleRate ;

EstimatedFilter = fliplr(cc) ;

% PREDICTING RESPONSE
EstimatedResponse = conv(Stim,EstimatedFilter) ;
EstimatedResponse = EstimatedResponse(MaxLag*sampleRate+1:end-MaxLag*sampleRate) ;

% figures

figure %model
subplot(3,1,1)
plot(time,Stim)
xlabel('time')
ylabel('Stim units')

subplot(3,1,2)
plot(time,Filter,'r')
xlabel('time')
ylabel('Filter units')

subplot(3,1,3)
plot(time,Response,'g')
xlabel('time')
ylabel('Response units')

figure % cc
plot(ccTime,cc)
xlabel('time')
ylabel('cc')

figure % filter estimate
plot(time,Filter) ;
hold on
plot(time,EstimatedFilter,'r--')
xlabel('time')
ylabel('Filter units')
legend('true', 'estimated')

figure % estimated response
plot(time,Response/max(Response))
hold on
plot(time,EstimatedResponse,'r--')
xlabel('time')
ylabel('response units')
legend('true', 'estimated')


