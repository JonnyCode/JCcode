% correlation and convolution simulation, This simulation will explore
% finding filters and prediciting responses in a purely linear system.
% It is similar to "FollowUpSim" but reduced

% JC 12/5/12

% STIMULUS

% stim and data collection parameters (parameters of experiment) 
sampleRate = 10000 ; % hz
MaxTime = 1 ; % sec

% time 
time = [1/sampleRate:1/sampleRate:MaxTime] ;

% stimulus - This is the stimulus you are probing the system with, I provided gaussian noise but you can replace this
% function with whatever you want to probe the system with and see how that helps or hurts your ability to detect the Filter

Stim = normrnd(0,1,1,length(time)) ; % gaussian noise

% MODEL

% model parameters (parameters of biology)
tpeak = 0.05 ; % sec (gaussian filtered 0 is centered, acausal, +X becomes more causal)
peakRise = .02 ; % sec

% Filter - This is the filter the we are pretending the biology implements (you can replace the gaussian with any function you want). 
FilterTime = time-(MaxTime/2) ;
Filter = exp(-((FilterTime-tpeak).^2)/(2*peakRise^2)) ; % gaussian

% linear response - We are pretending the system is perfectly linear at this point, so the response is just the convolution of the filter and stimulus
Response = conv(Stim,Filter,'same') ;

% ANALYSIS

% finding filter - As the expereimentor you only get the stimulus you provided and the response you got and you are trying to find the linear filter
cc= xcorr(Stim,Response,'coeff') ;
ccTime = ([1:length(cc)] - (length(cc)+1)/2)/sampleRate ;

if abs(min(cc))>abs(max(cc)) ;
    EstimatedFilter = fliplr(cc)/min(cc) ;
else
    EstimatedFilter = fliplr(cc)/max(cc) ;
end

% estimating response
EstimatedResponse = conv(Stim,EstimatedFilter,'same') ;

% figures

figure %model
subplot(3,1,1)
plot(time,Stim)
xlabel('time')
ylabel('Stim units')

subplot(3,1,2)
plot(FilterTime,Filter,'r')
xlabel('time')
ylabel('Filter units')

subplot(3,1,3)
plot(time,Response/max(Response),'g')
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

figure % estimated response
plot(time,Response/max(Response))
hold on
plot(time,EstimatedResponse/max(EstimatedResponse),'r--')
xlabel('time')
ylabel('response units')
legend('true', 'estimated')
