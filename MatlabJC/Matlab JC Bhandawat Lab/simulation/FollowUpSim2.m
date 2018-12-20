% correlation and convolution simulation, This simulation will explore
% finding filters and prediciting responses in a linear system with static nonlinearities.
% It is similar to "FollowUpSim" but has added nonlinearity

% JC 12/3/12

% STIMULUS

% stim and data collection parameters (parameters of experiment) 
sampleRate = 10000 ; % hz
MaxTime = 100 ; % sec

PulseTime = 2 ; % sec

StimType = 'noise' ; % 'Pulse' or nothing if you want gaussian noise
lpfCutOff = 5000; % hz , if you want to slow your probe stim down

% time 
time = [1/sampleRate:1/sampleRate:MaxTime] ;

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


% MODEL

% model parameters (parameters of biology)
tpeak = 0.05 ; % sec (gaussian filtered 0 is centered, acausal, +X becomes more causal)
peakRise = .02 ; % sec

% Filter - This is the filter the we are pretending the biology implements (you can replace the gaussian with any function you want). 
FilterTime = time-(MaxTime/2) ;
Filter = exp(-((FilterTime-tpeak).^2)/(2*peakRise^2)) ; % gaussian

% static nonlinearity - also implemented by biology (replace with any
% function you want, x axis will run from trough to peak of
% linear response)

%Nonlinearity = [-100:100] ; % linear
Nonlinearity = [-100:100].^2 ; % parabola
%Nonlinearity = cumsum(exp((-[-100:100].^2)/(30^2))) ; % aigmoid

% linear response - We are pretending the system is perfectly linear at this point, so the response is just the convolution of the filter and stimulus
LinearResponse = conv(Stim,Filter,'same') ;

% nonlinear response
InputX = ([0:length(Nonlinearity)-1]/(length(Nonlinearity)-1))*(max(LinearResponse)+abs(min(LinearResponse)))-abs(min(LinearResponse)) ; % scale x-axis
Response = interp1(InputX,Nonlinearity,LinearResponse,'linear','extrap') ; % implement nonlinearity


% ANALYSIS

% analysis parameters
MaxLag = 1 ; % sec, time length of filter you are finding
NumEstNLBins = 10 ; % number of bins to estimate nonlinearity

% finding filter - As the expereimentor you only get the stimulus you provided and the response you got and you are trying to find the linear filter
cc= xcorr(Stim,Response,MaxLag*sampleRate,'coeff') ;
ccTime = ([1:length(cc)] - (length(cc)+1)/2)/sampleRate ;

if abs(min(cc))>abs(max(cc)) ;
    EstimatedFilter = fliplr(cc)/min(cc) ;
else
    EstimatedFilter = fliplr(cc)/max(cc) ;
end

% estimating response
EstimatedLinearResponse = conv(Stim,EstimatedFilter) ;
EstimatedLinearResponse = EstimatedLinearResponse(MaxLag*sampleRate+1:end-MaxLag*sampleRate) ;

EstimatedNonlinearityX = [min(EstimatedLinearResponse):range(EstimatedLinearResponse)/NumEstNLBins:max(EstimatedLinearResponse)] ;
for a=1:length(EstimatedNonlinearityX)-1 ; % for each bin
    iEstimatedLinearResponse = find(EstimatedLinearResponse>=EstimatedNonlinearityX(a) & EstimatedLinearResponse<=EstimatedNonlinearityX(a+1)) ;
    EstimatedNonlinearity(a) = mean(Response(iEstimatedLinearResponse)) ;
end
EstimatedNonlinearityX = EstimatedNonlinearityX(1:end-1)+diff(EstimatedNonlinearityX)/2 ; % take center of bins

EstimatedResponse = interp1(EstimatedNonlinearityX,EstimatedNonlinearity,EstimatedLinearResponse,'linear','extrap') ;


% figures

figure %model
subplot(4,1,1)
plot(time,Stim)
xlabel('time')
ylabel('Stim units')

subplot(4,1,2)
plot(FilterTime,Filter,'r')
xlabel('time')
ylabel('Filter units')

subplot(4,1,3)
plot(Nonlinearity)
xlabel('input')
ylabel('output')

subplot(4,1,4)
plot(time,LinearResponse/max(LinearResponse),'g')
hold on
plot(time,Response/max(Response),'k')
xlabel('time')
ylabel('Response units')

figure % cc
subplot(2,1,1)
plot(time,Stim/max(Stim))
hold on
plot(time,Response/max(Response),'k')

subplot(2,1,2)
plot(ccTime,cc)
xlabel('time')
ylabel('cc')

figure % linear filter estimate
plot(FilterTime,Filter) ;
hold on
plot(ccTime,EstimatedFilter,'r--')
xlabel('time')
ylabel('Filter units')
legend('true', 'estimated')

figure % estimated nonlinearity
subplot(1,2,1)
plot(Nonlinearity)
xlabel('Input')
ylabel('Output')
title('true')

subplot(1,2,2)
plot(EstimatedLinearResponse,Response,'*')
hold on
plot(EstimatedNonlinearityX,EstimatedNonlinearity,'r')
xlabel('Linear Estimate')
ylabel('Actual Response')
legend('data','estimated nonlinearity')

figure % estimated response
plot(time,Response/max(Response))
hold on
plot(time,EstimatedResponse/max(EstimatedResponse),'r--')
xlabel('time')
ylabel('response units')
legend('true', 'estimated')


