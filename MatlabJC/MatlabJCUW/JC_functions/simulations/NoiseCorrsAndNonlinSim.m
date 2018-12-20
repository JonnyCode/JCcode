% simulation exploring advantage/disadvantages of nonlinearities in
% a circuit.  

% model: stimulus + common noise diverge into parrallel and othogonal
% information channels, nonlinearities in each channel are idealized to
% separate signal and noise, indepedant noise sources are added after
% linear nonlinear stage.

% questions: how does info retained by two channels depend on presence or
% absence of nonlinearity?  how does this change if channels are not
% orthogonal?

% Cloudy thinking - this model should well describe why noise correlations are stim dependant.
% nonlinearity will cut out common noise which is good bc noise is bad for info,
% but is bad for info bc common noise is better than independant noise.
% But in the case that addative noise is constant than nonlinearities
% should still help.

% parameters and logicals:
NumPnts = 1000000 ;

stimulus_std = 5 ;
stimulus_mean = 10 ;

noise_common_std = 5 ;
noise_common_mean = 5 ;

noise_ind_std = 1 ;
noise_ind_mean = 5 ;

OrthogonalizeFilters = 0 ;  % if you want orthogonal filters
AntiCorrFilters = 0 ; % if you want antiCorrelating filters
SameFilters = 1 ; % if you want channels to be the same
% all zeros will make filter2 the diriv of filter 1

% linear filters:

% filter 1 (only slightly more complicated than simple sum of gaussians)
filter1 = simFilter([1:NumPnts],50,30,10,100,50,1) ;  % simFilter(time,tpeak,peakRise,peakAmp,ttrough,troughDecay,troughAmp)
filter1 = filter1/sum(filter1) ;

% derivative of filter 1 (does this make channel correlations between 1 and 0?)
filter2 = [diff(filter1),0] ; 
filter2 = filter2/sum(filter2) ;

% grahm-schmit orthogonilization (channel signals seem to be independant only at Tau=0 but not otherwise)
if OrthogonalizeFilters == 1 ;
    filter2 = filter2 - filter1*((filter1*filter2')/(filter1*filter1')) ;
    filter2 = filter2/sum(filter2) ;
end

% if anticorr (does this make channel signals anticorrelated?)
if AntiCorrFilters == 1 ;
    filter2 = -filter1 ;
end

% if same filters (should makes signals perfectly correlated)
if SameFilters == 1 ;
    filter2 = filter1 ;
end

% stimulus and noise:
stimulus = normrnd(stimulus_mean,stimulus_std,1,NumPnts) ;


Noise_common = normrnd(noise_common_mean,noise_common_std,1,NumPnts) ;
Noise_1 = normrnd(noise_ind_mean,noise_ind_std,1,NumPnts) ;
Noise_2 = normrnd(noise_ind_mean,noise_ind_std,1,NumPnts) ;

stimulus = lowPassFilter(stimulus,10000,60) ; % (stim,samplerate,freqcuttoff)
Noise_common = lowPassFilter(Noise_common,10000,60) ;
Noise_1 = lowPassFilter(Noise_1,10000,60) ;
Noise_2 = lowPassFilter(Noise_2,10000,60) ;

% linear filter stimulus only and assess correlation:
Channel_1_sig = ifft(fft(stimulus).*fft(filter1)) ;
Channel_2_sig = ifft(fft(stimulus).*fft(filter2)) ;

cc_sig = xcov(Channel_1_sig,Channel_2_sig,'coef') ;

% linear filter noise only:
Channel_1_noise_common = ifft(fft(Noise_common).*fft(filter1)) ;
Channel_2_noise_common = ifft(fft(Noise_common).*fft(filter2)) ;

% linear filter stimulus and common noise (matters if added pre or post linear filter, so where is most interesting?):
Channel_1_plusNoise = ifft(fft(stimulus+Noise_common).*fft(filter1)) ;
Channel_2_plusNoise = ifft(fft(stimulus+Noise_common).*fft(filter2)) ;

% optimal nonlinearities for noise filtering by amplitude:
bins = [min([(noise_common_mean-2*noise_common_std),(stimulus_mean-2*stimulus_std)])...
    :max([(noise_common_mean+2*noise_common_std),(stimulus_mean+2*stimulus_std)])] ;
Channel1_sig_dist = hist(Channel_1_sig,bins) ;
Channel2_sig_dist = hist(Channel_2_sig,bins) ;
Channel1_noise_dist = hist(Channel_1_noise_common,bins) ;
Channel2_noise_dist = hist(Channel_2_noise_common,bins) ;

Channel1_sig_dist = Channel1_sig_dist/sum(Channel1_sig_dist) ; % prob dist
Channel2_sig_dist = Channel2_sig_dist/sum(Channel2_sig_dist) ;
Channel1_noise_dist = Channel1_noise_dist/sum(Channel1_noise_dist) ;
Channel2_noise_dist = Channel2_noise_dist/sum(Channel2_noise_dist) ;
 
% Nonlin1 = Channel1_dist./(Channel1_dist+Noise_dist) ;
% Nonlin2 = Channel2_dist./(Channel2_dist+Noise_dist) ;

Nonlin1 = zeros(1,length(bins)) ;
Nonlin2 = zeros(1,length(bins)) ;
Nonlin1(14:end) = 1 ;
Nonlin2(14:end) = 1 ;

% apply nonlinear weighting
Channel_1_postNL = nans(1,NumPnts) ;
Channel_2_postNL = nans(1,NumPnts) ;

for a=1:NumPnts ;
     i1 = find(bins==round(Channel_1_plusNoise(a))) ;
     
     Channel_1_postNL(a) = Channel_1_plusNoise(a)*Nonlin1(i1) ;
    
     i2 = find(bins==round(Channel_2_plusNoise(a))) ;
     Channel_2_postNL(a) = Channel_2_plusNoise(a)*Nonlin2(i2) ;
end

% add independent noise
Channel_1 = Channel_1_postNL + Noise_1 ;
Channel_2 = Channel_2_postNL + Noise_2 ;

% calculate information about stimulus from 2 channels (alone and in tandem)

% stim dependance of noise correlations
geoMean = sqrt(Channel_1.*Channel_2) ;
bins2 = min(geoMean):max(geoMean) ;
NoiseCorr = zeros(1,length(bins2)-1) ;

for a=1:length(bins2)-1 ;
    i1 = find(geoMean>bins2(a) & geoMean<=bins2(a+1)) ;

    if length(i1)>1 ;
        temp = corrcoef(Channel_1(i1)-Channel_1_sig(i1),Channel_2(i1)-Channel_2_sig(i1)) ;
        NoiseCorr(a) = temp (2,1) ;
    end
end


% figure
figure
subplot(2,1,1)
plot(bins,Channel1_sig_dist)
hold on
plot(bins,Channel1_noise_dist,'r')
plot(bins,Nonlin1,'k')

subplot(2,1,2)
plot(bins,Channel2_sig_dist)
hold on
plot(bins,Channel2_noise_dist,'r')
plot(bins,Nonlin2,'k')

figure
plot(bins2(1:length(NoiseCorr)),NoiseCorr)


