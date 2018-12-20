% this script was modified from "correlated_noise3 ". 


% SIMULATED CURRENT RECORDINGS
V = -30 ; % intermediate holding potential
rp_exc = 0 ; % exc rev potential
rp_inh = -60 ; %inh rev pot
err = 5 ;%mV  you held cell with error of reversal 

anticorrelated = 1 ; % if the noise added is simulating an anticorrelation in the conductances put a 1 here

meanStable_exc = 800 ; % mean of stable conductance sine wave
meanStable_inh = 800 ;

recFactor_exc = 0 ; % a negative will shift the sine wave into a rectification and positive will bring it out, zero means some of the noise will still be rectified
recFactor_inh = 0 ;


b=1 ;
for  sigma_common = 0:10:100 ; % std of common noise
    
a=1 ;
for hp = [-60,V,0] ; % for each holding potential
    
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
G_stable_exc = sin(.001*[1:lengthG])*meanStable_exc + meanStable_exc + recFactor_exc ; % stable exc g
G_stable_inh = sin(.001*[1:lengthG]+(45*pi))*meanStable_inh + meanStable_inh + recFactor_inh ; % stable inh g

% add noise to stable conductances
G_exc = G_stable_exc + G_noise2_exc ;
G_inh = G_stable_inh + G_noise2_inh ;

% rectify conductances
recPnts_exc = find(G_exc < 0) ; 
recPnts_inh = find(G_inh < 0) ;

G_exc(recPnts_exc) = 0 ;
G_inh(recPnts_inh) = 0 ;

% calculate recorded currents
I{a}(trial,1:lengthG) = G_exc*(hp-rp_exc) + G_inh*(hp-rp_inh) ;

end
a = a+1 ;
end

% ANALYZE SIMULATDED CURRENTS

% calculate standard deviation and variance of currents
var_Iexc(b) = mean(var(I{1},0,2)) ; 
var_Iint(b) = mean(var(I{2},0,2)) ; 
var_Iinh(b) = mean(var(I{3},0,2)) ; 

% calculate I residuals 
mean_Iexc{b} = mean(I{1}) ; 
mean_Iint{b} = mean(I{2}) ;
mean_Iinh{b} = mean(I{3}) ;

resid_Iexc{b} = I{1} - repmat(mean_Iexc{b},num_trials,1) ; % repmat replicates vector of mean
resid_Iint{b} = I{2} - repmat(mean_Iint{b},num_trials,1) ;
resid_Iinh{b} = I{3} - repmat(mean_Iinh{b},num_trials,1) ;

% variance of residual exc, inh, and exc+inh currents
var_ResIexc(b) = mean(var(resid_Iexc{b},0,2)) ;
var_ResIint(b) = mean(var(resid_Iint{b},0,2)) ;
var_ResIinh(b) = mean(var(resid_Iinh{b},0,2)) ;

% change exc and inh currents into conductances
G_fromI_exc{b}= I{1}./(rp_inh-rp_exc) ;
G_fromI_inh{b}= I{3}./(rp_exc-rp_inh) ;

% calculate G residuals 
mean_GfromI_exc{b} = mean(G_fromI_exc{b}) ; 
mean_GfromI_inh{b} = mean(G_fromI_inh{b}) ;

resid_GfromI_exc{b} = G_fromI_exc{b} - repmat(mean_GfromI_exc{b},num_trials,1) ; % repmat replicates vector of mean
resid_GfromI_inh{b} = G_fromI_inh{b} - repmat(mean_GfromI_inh{b},num_trials,1) ;

% variance of conductance residuals
var_ResGexc(b) = mean(var(resid_GfromI_exc{b},0,2)) ;
var_ResGinh(b) = mean(var(resid_GfromI_inh{b},0,2)) ;

% extract the factors that change Gexc and Ginh into Iint
factors = pinv([mean_GfromI_exc{b};mean_GfromI_inh{b}]')*mean_Iint{b}' ; % pinv is the pseudo inverse of the matrix because it is not a square matrix
ExcFact(b) = factors(1) ;
InhFact(b) = factors(2) ;

predicted_Iint{b} = ExcFact(b)*mean_GfromI_exc{b} + InhFact(b)*mean_GfromI_inh{b} ; 

% calculate the predicted variance of intermediate current based on the
% factors found above
Pred_varResIint(b) = (var_ResGexc(b)*(ExcFact(b)^2)) + (var_ResGinh(b)*(InhFact(b)^2)) ;


% Including power spectrum

% power spectrum of residual currents
[powerspec_xvalues_ResIexc{b}, mean_powerspec_ResIexc{b}] = PowerSpectrumFinder(resid_Iexc{b},10000) ; % power spectrum of exc residual currents
[powerspec_xvalues_ResIint{b}, mean_powerspec_ResIint{b}] = PowerSpectrumFinder(resid_Iint{b},10000) ; % power spectrum of int residual currents
[powerspec_xvalues_ResIinh{b}, mean_powerspec_ResIinh{b}] = PowerSpectrumFinder(resid_Iinh{b},10000) ; % power spectrum of inh residual currents

% power spectrum of residual conductances
[powerspec_xvalues_ResGexc{b}, mean_powerspec_ResGexc{b}] = PowerSpectrumFinder(resid_GfromI_exc{b},10000) ; % power spectrum of exc residual currents
[powerspec_xvalues_ResGinh{b}, mean_powerspec_ResGinh{b}] = PowerSpectrumFinder(resid_GfromI_inh{b},10000) ; % power spectrum of inh residual currents

% calculate the predicted power spectrum of the intermediate current as done above
pred_powerspecResIint{b} = (mean_powerspec_ResGexc{b}*(ExcFact(b)^2))+(mean_powerspec_ResGinh{b}*(InhFact(b)^2)) ;

PredvActual_powerResIint(b,:) = pred_powerspecResIint{b}./mean_powerspec_ResIint{b} ; % ratio of the predicted over the simulated
smoothed_powerspec_ratio(b,:) = SmoothPowerSpectrum(PredvActual_powerResIint(b,:), powerspec_xvalues_ResGexc{b}, 4, 1) ;% function to smooth powerspec for log scaling

b = b+1 ;
end

% ratio of predicted vs actual and confidence ratios for assumed error (amount you are off in voltage hp) 
PredvActual_varResIint = Pred_varResIint./var_ResIint ; % ratio of prediction vs actuall variance of intermediate current

assumed_drive = rp_exc-rp_inh ;
confidence_corr = ((assumed_drive^2)+(2*assumed_drive*err)+(2*(err^2)))/assumed_drive^2 ; % above this you are confident of declaring correlation between Gexc and Ginh
confidence_anticorr = ((assumed_drive^2)-(2*assumed_drive*err)+(2*(err^2)))/assumed_drive^2 ; % below this you are confident of declaring anticorrelation between Gexc and Ginh


% Figures
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
plot([0:10:100],Pred_varResIint) ;
xlabel('std common noise')
ylabel('variance of predicted residual int current')

subplot(2,2,3)
plot([0:10:100],var_ResIint)
xlabel('std common noise')
ylabel('variances of residual Iint')
hold on
plot([0:10:100],Pred_varResIint,'r')
legend('simulated', 'analytical')

subplot(2,2,4)
plot([0:10:100],PredvActual_varResIint)
ylabel('predicted/simulated int current var')
xlabel('std common noise')
hold on
plot([0:10:100],ones(1,length([0:10:100]))*confidence_corr,'k')
plot([0:10:100],ones(1,length([0:10:100]))*confidence_anticorr,'k')

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
plot(powerspec_xvalues_ResIint{1},pred_powerspecResIint{1}(1:length(powerspec_xvalues_ResIint{1})),'r')
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
plot(powerspec_xvalues_ResIint{2},pred_powerspecResIint{2}(1:length(powerspec_xvalues_ResIint{2})),'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz)')
ylabel('power Residual Iint')
title('common noise std increases in each subplot')

subplot(2,2,2)
plot(powerspec_xvalues_ResIint{4},mean_powerspec_ResIint{4}(1:length(powerspec_xvalues_ResIint{4})))
hold on,
plot(powerspec_xvalues_ResIint{4},pred_powerspecResIint{4}(1:length(powerspec_xvalues_ResIint{4})),'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz)')
ylabel('power Residual Iint')

subplot(2,2,3)
plot(powerspec_xvalues_ResIint{8},mean_powerspec_ResIint{8}(1:length(powerspec_xvalues_ResIint{8})))
hold on,
plot(powerspec_xvalues_ResIint{8},pred_powerspecResIint{8}(1:length(powerspec_xvalues_ResIint{8})),'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz)')
ylabel('power Residual Iint')

subplot(2,2,4)
plot(powerspec_xvalues_ResIint{10},mean_powerspec_ResIint{10}(1:length(powerspec_xvalues_ResIint{10})))
hold on,
plot(powerspec_xvalues_ResIint{10},pred_powerspecResIint{10}(1:length(powerspec_xvalues_ResIint{10})),'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz)')
ylabel('power Residual Iint')


figure,
plot(powerspec_xvalues_ResIint{10},PredvActual_powerResIint(:,1:length(powerspec_xvalues_ResIint{10})))
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

