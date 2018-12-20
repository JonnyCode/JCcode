% This function is modified from "precsion model" and will use a set of
% gaussian conductances and their PSTH to determine the type of
% conductances that leads to the most sterotyped output. 

% This function will implement a basic integrate and fire model for a set of
% noisy conductances.  Hopefully, this model will help answer the
% question: What synaptic inputs will permit high spike precision regardless
% of intrinsic mechanisms?

%NOTE- sections between lines are helpful when running function with only a
%single set of conducatance parameters.

% J. Cafaro 2/4/08


% clear all
% close all


% set time (ms)  
tmax = 1000 ;           % max of time
h= .1 ;                 % time step 
time = h:h:tmax ;       % time points assessed
tpoints=length(time) ;  % number of time points

TrialNumber = 2 ;      % number of trial to produce spikes

% set conductance noise parameters 
noise_max = 10 ;%ns   % the point at which noise stops increasing as function of g amp
noise_min = 0 ;%ns      
noise_start = 0 ;%ns   % the amp of the conductance at which the noise begins growing 
noise_stop = 10 ;%ns    % the amp of the conductance at which the noise maxes out (sigmoid points)  

freqcutoffNoise_exc = 60 ; % Hz high frequency cutoff (highest frequency allowed)
freqcutoffNoise_inh = 60 ; % Hz
freqcutoffNoise_common = 60 ; %Hz

muNoise_inh = 0 ;%nS       % mean of noise for inhibitory conductance
muNoise_exc = 0 ;          % mean of noise for excitatory conductance

correlation =1 ; % 1=correlated or -1=anticorrelated common noise
correlation_lag = 3/h ;%ms  inhibitory lag of correlation

% set intrinsic cell properties
E_inh = -80 ;%mV     % reversal potential for inhibitory current
E_exc = 0 ;          % reversal potential for excitatory current
E_leak = -90 ;       % reversal potential of leak conductance (this can make cell spike spontaneously)

C = 0.08 ; %nF       % Set capacitance of cell
G_leak = 25 ; %nS    % Set leak coductance    
V_rest = -60 ;       % initial potential of cell (resting potential mV) 
V_thr = -40 ;        % spike threshold
  
abs_ref = 2/h ;      % absolute refractory period (2ms)

% set stable excitatory and inhibitory white noise conductances parameters
muStable_inh = 2 ; % the mean of the inhibitory conductance
muStable_exc = 2 ; % the mean of the excitatory conductance

sigmaStable_inh = 100 ; % the standard deviation
sigmaStable_exc = 100 ;

freqcutoffStable_exc = 40 ; % Hz high frequency cutoff (highest frequency allowed)
freqcutoffStable_inh = 40 ; % Hz

% create stable conductances...
preGstable_exc = normrnd(muStable_inh,sigmaStable_inh,1,tpoints) ; % make a random exc g
preGstable_inh = zeros(1,length(preGstable_exc)) ; %  make a zeros block for inh g
preGstable_inh(correlation_lag+1:end) = preGstable_exc(1:end-correlation_lag) * (sigmaStable_inh/sigmaStable_exc)*correlation ; % same as gexc but a different amplitude, sign, and lag

% filter them, and rectify them
freqcutoffStable_exc = freqcutoffStable_exc/((1000/h)/tpoints) ; 
GstableExc_fft = fft(preGstable_exc) ; % take fft of exc signal
GstableExc_fft(:,1+freqcutoffStable_exc:length(GstableExc_fft)-freqcutoffStable_exc) = 0 ; % cut out high frequencies in first and second half of fft
prerecGstable_exc = real(ifft(GstableExc_fft)) ; % inverse fft
negpnts = find(prerecGstable_exc<0); %find the indices of the negative points
Gstable_exc = prerecGstable_exc ;
Gstable_exc(negpnts) = 0 ;

freqcutoffStable_inh = freqcutoffStable_inh/((1000/h)/tpoints) ; %gets freq cuttoff into apropriate units
GstableInh_fft = fft(preGstable_inh) ; % take fft of exc signal
GstableInh_fft(:,1+freqcutoffStable_inh:length(GstableInh_fft)-freqcutoffStable_inh) = 0 ; % cut out high frequencies in first and second half of fft
prerecGstable_inh = real(ifft(GstableInh_fft)) ; % inverse fft
negpnts = find(prerecGstable_inh<0); %find the indices of the negative points
Gstable_inh = prerecGstable_inh ;
Gstable_inh(negpnts) = 0 ;

% create noise for conductance
sigmaNoise_inh = (smf(Gstable_inh,[noise_start,noise_stop])*(noise_max-noise_min))+noise_min ;  % standard deviation of noise for inhibitory conductance
sigmaNoise_exc = (smf(Gstable_exc,[noise_start,noise_stop])*(noise_max-noise_min))+noise_min ;  % these values derived from sigmoidal relationship between g-amp and noise amp
sigmaNoise_common = (sigmaNoise_inh+sigmaNoise_exc)./2 ; % the common noise will be the mean of the two other sources (think about this more)

sigmaNoise_inh = sigmaNoise_inh*(1-FracCommon_noise) ; % make the standard deviation reflect the stdev that will be added for the common noise (this is not wierd and should be better considered)
sigmaNoise_exc = sigmaNoise_exc*(1-FracCommon_noise) ; % this kind of code will mean the the variance of the noise will shrink as more of the noise is common noise
sigmaNoise_common = sigmaNoise_common*FracCommon_noise ;

Spikes = zeros(TrialNumber,tpoints) ; % prep spike matrix

for trial = 1:TrialNumber ; % for each trial

% reset variables to go for trial
V_trace(trial,:) = time ; 
V_trace(trial,1) = V_rest ;
I_syn = time ;
I_syn(1) = 0;

ref = 0 ;    
    
% set noise for this trial    
preGnoise_common = normrnd(0,sigmaNoise_common,1,tpoints) ;
freqcutoffNoise2_common = freqcutoffNoise_common/((1000/h)/tpoints) ; 
GnoiseCom_fft = fft(preGnoise_common) ; % take fft of exc signal
GnoiseCom_fft(:,1+freqcutoffNoise2_common:length(GnoiseCom_fft)-freqcutoffNoise2_common) = 0 ; % cut out high frequencies in first and second half of fft
prerecGnoise_common = real(ifft(GnoiseCom_fft)) ; % inverse fft

preGnoise_inh = normrnd(0,sigmaNoise_inh) ;
freqcutoffNoise2_inh = freqcutoffNoise_inh/((1000/h)/tpoints) ; 
GnoiseInh_fft = fft(preGnoise_inh) ; % take fft of exc signal
GnoiseInh_fft(:,1+freqcutoffNoise2_inh:length(GnoiseInh_fft)-freqcutoffNoise2_inh) = 0 ; % cut out high frequencies in first and second half of fft
prerecGnoise_inh = real(ifft(GnoiseInh_fft)) ; % inverse fft
Gnoise_inh(trial,:) = prerecGnoise_inh + prerecGnoise_common ;

preGnoise_exc = normrnd(0,sigmaNoise_exc) ;
freqcutoffNoise2_exc = freqcutoffNoise_exc/((1000/h)/tpoints) ; 
GnoiseExc_fft = fft(preGnoise_exc) ; % take fft of exc signal
GnoiseExc_fft(:,1+freqcutoffNoise2_exc:length(GnoiseExc_fft)-freqcutoffNoise2_exc) = 0 ; % cut out high frequencies in first and second half of fft
prerecGnoise_exc = real(ifft(GnoiseExc_fft)) ; % inverse fft
Gnoise_exc(trial,:) = prerecGnoise_exc + correlation*prerecGnoise_common ;

% add noise to the stable conductances and rectify
G_exc(trial,:) = Gstable_exc + Gnoise_exc(trial,:) ;
negpnts = find(G_exc(trial,:)<0) ;
G_exc(trial,negpnts) = 0 ;

G_inh(trial,:) = Gstable_inh + Gnoise_inh(trial,:) ;
negpnts = find(G_inh(trial,:)<0) ;
G_inh(trial,negpnts) = 0 ;

for t= 2 : tpoints

% LIF model function
I_syn(t) = G_exc(trial,t)*(V_trace(trial,t-1)-E_exc) + G_inh(trial,t)*(V_trace(trial,t-1)-E_inh) + G_leak*(V_trace(trial,t-1)-E_leak) ;

if ref == 0 ; % if not within refractory period
V_trace(trial, t) = V_trace(trial, t-1) + (h/1000)*(-I_syn(t)/C) ; %h needs to be in units of seconds from ms

else
V_trace(trial, t) = V_rest ;
ref = ref - 1 ;

end

if V_trace(trial,t)>V_thr
V_trace(trial,t) = 50 ;
ref = abs_ref ;
end

end % end t loop

% Note times of spikes
Spikepnts = find(V_trace(trial,:) == 50) ;
Spikes(trial,Spikepnts) = 1 ;

figure(1)
plot(V_trace(trial,:))
hold on
plot((Spikes(trial,:)*40)+V_thr,'r')

figure
subplot(4,1,1)
plot(Gnoise_inh(trial,:),'r')
subplot(4,1,2)
plot(Gnoise_exc(trial,:))
subplot(4,1,3)
plot(G_exc(trial,:))
hold on
plot(G_inh(trial,:),'r')
plot(Gstable_exc,'b-.')
plot(Gstable_inh,'r-.')
subplot(4,1,4)
plot(V_trace(trial,:))
hold on
plot((Spikes(trial,:)*40)+V_thr,'r')

end % this ends the trial loop

% assess spike precision
AllSpikepnts = find(sum(Spikes)>0) ; % find all the pnts where a spike occured on any trial

for a= 1:length(find(AllSpikepnts<tpoints-100 & AllSpikepnts) ; % for each spike that occured in all trials before the 100th to last point (so that we can still see the following G ....
NonperciseG_exc(a,:) = Gstable_exc(AllSpikepnts(a)-100:AllSpikepnts(a)+100) ;
NonperciseG_inh(a,:) = Gstable_inh(AllSpikepnts(a)-100:AllSpikepnts(a)+100) ;
end
clear a
meanNonpreciseG_exc = mean(NonperciseG_exc) ;
meanNonpreciseG_inh = mean(NonperciseG_inh) ;


PSTH = smooth(sum(Spikes),100) ; % smooths the PSTH by 10 time points
PSTH_peakpnt = find(PSTH==max(PSTH(100:end))) ; % find the most percise spikes times

if length(PSTH_peakpnt)>1            
    diff_PSTH_peakpnt = find(diff(PSTH_peakpnt)>1); % find all the non consecuative PSTH peaks
    if ~isempty(diff_PSTH_peakpnt) % if more than one highest peak point exists ...
        PSTH_peaks = PSTH_peakpnt([1,(diff_PSTH_peakpnt + 1)']) ; % get the PSTH_peakpnt that are the start of each PSTH peak 
    else % if only one highest peak point exists
        PSTH_peaks = PSTH_peakpnt(1) ;
    end
else PSTH_peaks = PSTH_peakpnt ;
end

for a = 1:length(PSTH_peaks) ; % for each "percise spike" ...
perciseG_exc(a,:) = Gstable_exc(PSTH_peaks(a)-100:PSTH_peaks(a)+100) ; % note the stable exc wave form that caused the spike 
perciseG_inh(a,:) = Gstable_inh(PSTH_peaks(a)-100:PSTH_peaks(a)+100) ;
diff_perciseG_exc(a,:) = perciseG_exc(a,:) - meanNonpreciseG_exc ; % what is the difference between the G that caused the most percise spikes and the mean G that caused all spikes
diff_perciseG_inh(a,:) = perciseG_inh(a,:) - meanNonpreciseG_inh ;
end
clear a




%__________________________________________________________________________
% plot conductances
figure
plot(Gstable_inh,'r')
hold on
plot(Gstable_exc)
title('base conductances')
plot(PSTH*10,'g')
legend('Ginh','Gexc','PSTH*10')

figure
plot([1:size(perciseG_exc,2)],perciseG_exc)
hold on
plot([1:size(perciseG_inh,2)],perciseG_inh,'--')
title('percise G')

figure
plot(meanNonpreciseG_exc)
hold on
plot(meanNonpreciseG_inh,'r')
title('nonpercise G')

figure
plot([1:size(diff_perciseG_exc,2)],diff_perciseG_exc)
hold on
plot([1:size(diff_perciseG_inh,2)],diff_perciseG_inh,'--')
title('difference of percise and non')

%________________________________________________________________________

end % this ends the function