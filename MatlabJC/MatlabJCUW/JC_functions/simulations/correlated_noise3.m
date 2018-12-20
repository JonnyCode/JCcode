% this script was modified from "correlated_noise2 ". 


% SIMULATED CURRENT RECORDINGS
V = -20 ; % intermediate holding potential
anticorrelated = 0 ; % if the noise added is simulating an anticorrelation in the conductances put a 1 here

b=1 ;
for  sigma_common = 0:10:100 ; % std of common noise
    
a=1 ;
for hp = [-60,V,20] ; % for each holding potential
    
num_trials = 10 ; % the number of trials    
    
for trial = 1:num_trials ;

% construct gaussian conductance noise
mu = 0 ;                  % mean
sigma_exc = 50 ;             % std
sigma_inh = 50 ;             % std
lengthG = 20000 ;         % length of conductance     
filter = 1 ;               % filter (1 negates filter)
G_noise_exc = conv(normrnd(mu,sigma_exc,1,lengthG),filter) ; % convolved w/ filter
G_noise_inh = conv(normrnd(mu,sigma_inh,1,lengthG),filter) ;
G_common_noise = conv(normrnd(mu,sigma_common,1,lengthG),filter) ;

% low pass filter G noise
freqcutoff_exc =  10000 ; %Hz max frequency allowed
freqcutoff_inh = 10000 ;
freqcutoff_common = 10000 ;

freqcutoff_exc = freqcutoff_exc/(10000/lengthG) ; % this adjusts the freq cutoff for the length and assumes a sample rate of 10kHz
freqcutoff_inh = freqcutoff_inh/(10000/lengthG) ;
freqcutoff_common = freqcutoff_common/(10000/lengthG) ;

GnoiseLP_fft_exc = fft(G_noise_exc) ; % take fft of exc signal
GnoiseLP_fft_exc(:,1+freqcutoff_exc:length(GnoiseLP_fft_exc)-freqcutoff_exc) = 0 ; % cut out high frequencies in first and second half of fft
GnoiseLP_exc = real(ifft(GnoiseLP_fft_exc)); % inverse fft

GnoiseLP_fft_inh = fft(G_noise_inh) ; % take fft of inh signal
GnoiseLP_fft_inh(:,1+freqcutoff_inh:length(GnoiseLP_fft_inh)-freqcutoff_inh) = 0 ; % cut out high frequencies in first and second half of fft
GnoiseLP_inh = real(ifft(GnoiseLP_fft_inh)); % inverse fft

GnoiseLP_fft_com = fft(G_common_noise) ; % take fft of common signal
GnoiseLP_fft_com(:,1+freqcutoff_common:length(GnoiseLP_fft_com)-freqcutoff_common) = 0 ; % cut out high frequencies in first and second half of fft
GnoiseLP_com = real(ifft(GnoiseLP_fft_com)); % inverse fft

% add common noise
if anticorrelated == 1 ; % if the noise added is simulating an anticorrelation in the conductances
    G_noise2_exc = GnoiseLP_exc + GnoiseLP_com ; % add the low pass filtered exc noise to the lp filtered common noise 
    G_noise2_inh = GnoiseLP_inh - GnoiseLP_com ; % same as above but with inh
else                    % otherwise assume the noise is correlated
    G_noise2_exc = GnoiseLP_exc + GnoiseLP_com ; % add the low pass filtered exc noise to the lp filtered common noise 
    G_noise2_inh = GnoiseLP_inh + GnoiseLP_com ; % same as above but with inh
end


% construct stable conductance
G_stable_exc = sin(.001*[1:lengthG])*800+800 ; % stable exc g
G_stable_inh = sin(.001*[1:lengthG]+(45*pi))*800+800 ; % stable inh g

% add noise to stable conductances
G_exc = G_stable_exc + G_noise2_exc ;
G_inh = G_stable_inh + G_noise2_inh ;

% calculate recorded currents
rp_exc = 20 ; % exc rev potential
rp_inh = -60 ; %inh rev pot
I{a}(trial,1:lengthG) = G_exc*(hp-rp_exc) + G_inh*(hp-rp_inh) ;

end
a = a+1 ;
end

% ANALYZE SIMULATDED CURRENTS

% calculate standard deviation and variance of currents
var_Iexc(b) = var(I{1}(trial,1:lengthG),[],2) ; 
var_Iint(b) = var(I{2}(trial,1:lengthG),[],2) ; 
var_Iinh(b) = var(I{3}(trial,1:lengthG),[],2) ; 

% calculate I residuals 
mean_Iexc = repmat(mean(I{1}),num_trials,1) ; % repmat replicates vector of mean
mean_Iint = repmat(mean(I{2}),num_trials,1) ;
mean_Iinh = repmat(mean(I{3}),num_trials,1) ;

resid_Iexc = I{1} - mean_Iexc ;
resid_Iint = I{2} - mean_Iint ;
resid_Iinh = I{3} - mean_Iinh ;

% variance of exc, inh, and exc+inh currents
var_ResIexc(b) = mean(var(resid_Iexc,0,2)) ;
var_ResIint(b) = mean(var(resid_Iint,0,2)) ;
var_ResIinh(b) = mean(var(resid_Iinh,0,2)) ;

% change exc and inh currents into conductances
G_fromI_exc= I{1}./(rp_inh-rp_exc) ;
G_fromI_inh= I{3}./(rp_exc-rp_inh) ;

% calculate G residuals 
mean_GfromI_exc = repmat(mean(G_fromI_exc),num_trials,1) ; % repmat replicates vector of mean
mean_GfromI_inh = repmat(mean(G_fromI_inh),num_trials,1) ;

resid_GfromI_exc = G_fromI_exc - mean_GfromI_exc ;
resid_GfromI_inh = G_fromI_inh - mean_GfromI_inh ;

% variance of conductance residuals
var_ResGexc(b) = mean(var(resid_GfromI_exc,0,2)) ;
var_ResGinh(b) = mean(var(resid_GfromI_inh,0,2)) ;

% calculate the predicted variance of intermediate current based on pseudo analytical equation
Pred_varIint(b) = (var_ResGexc(b)*((V-rp_exc)^2)) + (var_ResGinh(b)*((V-rp_inh)^2)) ;


% Including power spectrum

% power spectrum of residual currents
[powerspec_xvalues_ResIexc{b}, mean_powerspec_ResIexc{b}] = PowerSpectrumFinder(resid_Iexc,10000) ; % power spectrum of exc residual currents
[powerspec_xvalues_ResIint{b}, mean_powerspec_ResIint{b}] = PowerSpectrumFinder(resid_Iint,10000) ; % power spectrum of int residual currents
[powerspec_xvalues_ResIinh{b}, mean_powerspec_ResIinh{b}] = PowerSpectrumFinder(resid_Iinh,10000) ; % power spectrum of inh residual currents

% power spectrum of residual conductances
[powerspec_xvalues_ResGexc{b}, mean_powerspec_ResGexc{b}] = PowerSpectrumFinder(resid_GfromI_exc,10000) ; % power spectrum of exc residual currents
[powerspec_xvalues_ResGinh{b}, mean_powerspec_ResGinh{b}] = PowerSpectrumFinder(resid_GfromI_inh,10000) ; % power spectrum of inh residual currents

% calculate the predicted power spectrum of the intermediate current as done above
pred_powerspecIint{b} = (mean_powerspec_ResGexc{b}*((V-rp_exc)^2))+(mean_powerspec_ResGinh{b}*((V-rp_inh)^2)) ;

PredvActual_powerIint(b,:) = pred_powerspecIint{b}./mean_powerspec_ResIint{b} ; % ratio of the predicted over the simulated
smoothed_powerspec_ratio(b,:) = SmoothPowerSpectrum(PredvActual_powerIint(b,:), powerspec_xvalues_ResGexc{b}, 4, 1) ;% function to smooth powerspec for log scaling

b = b+1 ;
end

PredvActual_varIint = Pred_varIint./var_ResIint ; % ratio of prediction vs actuall variance of intermediate current


%Figures
figure, plot([0:10:100],var_ResIexc) ;
xlabel('std common noise')
ylabel('variance of residual exc current')


figure, plot([0:10:100],var_ResIinh) ;
xlabel('std common noise')
ylabel('variance of residual inh current')


figure, subplot(2,2,1),
plot([0:10:100],var_ResIint) ;
xlabel('std common noise')
ylabel('variance of actual residual int current')

subplot(2,2,2)
plot([0:10:100],Pred_varIint) ;
xlabel('std common noise')
ylabel('variance of predicted residual int current')

subplot(2,2,3)
plot([0:10:100],var_ResIint)
xlabel('std common noise')
ylabel('variances of residual Iint')
hold on
plot([0:10:100],Pred_varIint,'r')
legend('simulated', 'analytical')

subplot(2,2,4)
plot([0:10:100],PredvActual_varIint)
ylabel('predicted/simulated int current var')
xlabel('std common noise')


figure
subplot(2,2,1)
plot(powerspec_xvalues_ResGexc{1},mean_powerspec_ResGexc{1}(1:length(powerspec_xvalues_ResGexc{1})))
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz)')
ylabel('power Residual Gexc')

subplot(2,2,2)
plot(powerspec_xvalues_ResGinh{1},mean_powerspec_ResGinh{1}(1:length(powerspec_xvalues_ResGinh{1})))
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz)')
ylabel('power Residual Ginh')

subplot(2,2,3)
plot(powerspec_xvalues_ResIint{1},mean_powerspec_ResIint{1}(1:length(powerspec_xvalues_ResIint{1})))
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz)')
ylabel('power Residual Iint')

subplot(2,2,4)
plot(powerspec_xvalues_ResIint{1},mean_powerspec_ResIint{1}(1:length(powerspec_xvalues_ResIint{1})))
hold on,
plot(powerspec_xvalues_ResIint{1},pred_powerspecIint{1}(1:length(powerspec_xvalues_ResIint{1})),'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz)')
ylabel('power Residual Iint')
legend('actual','predicted')
title('no common noise added')

figure
subplot(2,2,1)
plot(powerspec_xvalues_ResIint{2},mean_powerspec_ResIint{2}(1:length(powerspec_xvalues_ResIint{2})))
hold on,
plot(powerspec_xvalues_ResIint{2},pred_powerspecIint{2}(1:length(powerspec_xvalues_ResIint{2})),'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz)')
ylabel('power Residual Iint')
title('common noise std increases in each subplot')

subplot(2,2,2)
plot(powerspec_xvalues_ResIint{4},mean_powerspec_ResIint{4}(1:length(powerspec_xvalues_ResIint{4})))
hold on,
plot(powerspec_xvalues_ResIint{4},pred_powerspecIint{4}(1:length(powerspec_xvalues_ResIint{4})),'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz)')
ylabel('power Residual Iint')

subplot(2,2,3)
plot(powerspec_xvalues_ResIint{8},mean_powerspec_ResIint{8}(1:length(powerspec_xvalues_ResIint{8})))
hold on,
plot(powerspec_xvalues_ResIint{8},pred_powerspecIint{8}(1:length(powerspec_xvalues_ResIint{8})),'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz)')
ylabel('power Residual Iint')

subplot(2,2,4)
plot(powerspec_xvalues_ResIint{10},mean_powerspec_ResIint{10}(1:length(powerspec_xvalues_ResIint{10})))
hold on,
plot(powerspec_xvalues_ResIint{10},pred_powerspecIint{10}(1:length(powerspec_xvalues_ResIint{10})),'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz)')
ylabel('power Residual Iint')


figure,
plot(powerspec_xvalues_ResIint{10},PredvActual_powerIint(:,1:length(powerspec_xvalues_ResIint{10})))
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz)')
ylabel('ratio of predicted power Residual Iint/ actual')
legend('0','10','20','30','40','50','60','70','80','90','100')

figure,
plot(smoothed_powerspec_ratio(1,:).Freq,...
     [smoothed_powerspec_ratio(1,:).PowerSpec;smoothed_powerspec_ratio(2,:).PowerSpec;...
     smoothed_powerspec_ratio(3,:).PowerSpec;smoothed_powerspec_ratio(4,:).PowerSpec;...
     smoothed_powerspec_ratio(5,:).PowerSpec;smoothed_powerspec_ratio(6,:).PowerSpec;...
     smoothed_powerspec_ratio(7,:).PowerSpec;smoothed_powerspec_ratio(8,:).PowerSpec;...
     smoothed_powerspec_ratio(9,:).PowerSpec;smoothed_powerspec_ratio(10,:).PowerSpec;...
     smoothed_powerspec_ratio(11,:).PowerSpec])
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
xlabel('frequency (Hz)')
ylabel('ratio')
title('ratio of predicted power Residual Iint/ actual for each stdev of common noise')
legend('0','10','20','30','40','50','60','70','80','90','100')

