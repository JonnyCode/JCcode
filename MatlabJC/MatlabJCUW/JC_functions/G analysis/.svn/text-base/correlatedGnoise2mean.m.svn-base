function [Idata] = correlatedGnoise2mean(CellInfo_str,epochCond_num_str,exc_epochs_str, inh_epochs_str, int_epochs_str, firstsectiontv_str,lastsectiontv_str,firstsectionmean_str,lastsectionmean_str) ;

% This function will analyze conductances recordings and create figures in the same way as correlatedGnoise2 
% but will do this for signals with a mean light that follows a time varying signal so that it can use the Exc
% and Inh factors to change excitation and inhibition into intermediate current.

% Input = CellInfo file name in Index, epochCondtion number, epoch numbers for exc, inh, and
% intermediate currents, and the range of the sample points to be analyzed.
% ENTER ALL INPUT AS STRINGS:  eg ~ exc_epochs_str = '[5:10, 30:35]'
% Output = some numbers that may be helpful


% JC 7/27/07
% JC 1/08

% % UNCOMMMENT IF RUNNING AS SCRIPT
% % epochCondition
% CellInfo_str = '080405c10_a' ;
% epochCond_num = 1 ; % which epochcondition 
% 
% % Epochs recorded at exc, inh, and intermediate potentials
% exc_epochs = [531:540] ;
% inh_epochs = [464:473] ;
% int_epochs = [577:586] ;
% 
% % section of epoch to be analyzed (by sample points)
% firstsection = 2000;       %2000 or 64000
% lastsection = 75000;       %61000 or 121000


% FORMAT DATA

% load CellInfo file
cd ~/data_analysis/Index ;
load (CellInfo_str) ;  % parathesies allow it to load that varaible name

% Get appropriate cell structure format
CellInfo = LoadSCIData(CellInfo,1) ;  % for two amps used

EpochCondition = LoadAndSmoothEpochCondition(CellInfo,5000,1) ; 

% change strings into numbers
epochCond_num = str2num(epochCond_num_str) ;
exc_epochs = str2num(exc_epochs_str) ;
inh_epochs = str2num(inh_epochs_str) ;
int_epochs = str2num(int_epochs_str) ;
firstsectiontv = str2num(firstsectiontv_str) ;
lastsectiontv = str2num(lastsectiontv_str) ;
firstsectionmean = str2num(firstsectionmean_str) ;
lastsectionmean = str2num(lastsectionmean_str) ;


% PARAMETERS
% driving forces and holding potentials
E_Ampa = 0 ; % reversal potential for ampa channels
E_Gaba = -60 ; % reversal potential for gaba channel
exc_drive = E_Gaba - E_Ampa ; % E_holding - E_ampa  
inh_drive = E_Ampa - E_Gaba ; % E_holding - E_Gaba
int_hold = E_Ampa + E_Gaba /2 ;

% sample rate
samplerate = 10000 ; %Hz

% high frequency cutoff (highest frequency allowed)
freqcutoff_exc = 500 ; % Hz 
freqcutoff_inh = 500 ; % Hz
freqcutoff_int = 500 ; % Hz

% miscillaneous
amp_bins = [0:100] ; % bins for the conductance amplitude histogram

% FROM PARAMETERS
% change epoch numbers into epochcondition.epochdata.data row indices
for a = 1:length(exc_epochs) ;
    exc_data(a) = find(EpochCondition(epochCond_num).EpochNumbers == exc_epochs(a)) ;
end
for b = 1:length(inh_epochs) ;
    inh_data(b) = find(EpochCondition(epochCond_num).EpochNumbers == inh_epochs(b)) ;
end
for c = 1:length(int_epochs) ;
    int_data(c) = find(EpochCondition(epochCond_num).EpochNumbers == int_epochs(c)) ;
end
clear a b c 

% number of traces to be analyzed
numexc = length(exc_data) ;
numinh = length(inh_data) ;
numint = length(int_data) ; 
numtraces = numexc + numinh + numint ; 


% GET EXC AND INH FACTORS FROM TIME VARYING DATA SECTIONS
% separate data by holding potential and keep only section to be analyzed
Idata{1} = EpochCondition(epochCond_num).EpochData.Data(exc_data,firstsectiontv:lastsectiontv) ;  
Idata{2} = EpochCondition(epochCond_num).EpochData.Data(inh_data,firstsectiontv:lastsectiontv) ;
Idata{3} = EpochCondition(epochCond_num).EpochData.Data(int_data,firstsectiontv:lastsectiontv) ;

% find offset for each exc and inhibitory trace to assume isolation
Idata_offsets{1} = max(EpochCondition(epochCond_num).EpochData.Data(exc_data,:),[],2) ;  
Idata_offsets{2} = min(EpochCondition(epochCond_num).EpochData.Data(inh_data,:),[],2) ;

% adjust exc and inh current data by their offsets  
Idata{1} = Idata{1} - repmat(Idata_offsets{1},1,length(Idata{1})) ;
Idata{2} = Idata{2} - repmat(Idata_offsets{2},1,length(Idata{2})) ;

% low pass filter currents to remove noise from instrumentation
freqcutoff_exc = freqcutoff_exc/(10000/length(Idata{1})) ; % this adjusts the freq cutoff for the length and assumes a sample rate of 10kHz
freqcutoff_inh = freqcutoff_inh/(10000/length(Idata{2})) ;
freqcutoff_int = freqcutoff_int/(10000/length(Idata{3})) ;

Iexc_fft = fft(Idata{1},[],2) ; % take fft of exc signal
Iexc_fft(:,1+freqcutoff_exc:length(Iexc_fft)-freqcutoff_exc) = 0 ; % cut out high frequencies in first and second half of fft
Idata{1} = real(ifft(Iexc_fft,[],2)) ; % inverse fft

Iinh_fft = fft(Idata{2},[],2) ; % take fft of inh signal
Iinh_fft(:,1+freqcutoff_inh:length(Iinh_fft)-freqcutoff_inh) = 0 ; % cut out high frequencies in first and second half of fft
Idata{2} = real(ifft(Iinh_fft,[],2)) ; % inverse fft

Iint_fft = fft(Idata{3},[],2) ; % take fft of inh signal
Iint_fft(:,1+freqcutoff_int:length(Iint_fft)-freqcutoff_int) = 0 ; % cut out high frequencies in first and second half of fft
Idata{3} = real(ifft(Iint_fft,[],2)) ; % inverse fft

% CURRENTS
% calculate mean current at each holding potential
mean_Iexc = mean(Idata{1}) ;
mean_Iinh = mean(Idata{2}) ;
mean_Iint = mean(Idata{3}) ;

% CONDUCTANCES 
% calculate conductances
warning off last    % prevents NaNs warning about case sensitive bs
G_exc = NaNs(size(Idata{1}));   % prepares matrix
warning off last
G_inh = NaNs(size(Idata{2}));

for a = 1 : numexc ;  % for each exc current recorded
    G_exc(a,:) = Idata{1}(a,:)./exc_drive ; % calculate the exc conductance
end

for b = 1 : numinh ;  % for each inh current recorded
    G_inh(b,:) = Idata{2}(b,:)./inh_drive ; % calculate the exc conductance
end
clear a b

mean_Gexc = mean(G_exc) ;
mean_Ginh = mean(G_inh) ;

% calculate an intermediate current by finding two factors excFact and
% inhFact such that mean(excFact*Iexc + inhFact*Iinh - Iint)^2 is minimized  

factors = pinv([mean_Gexc;mean_Ginh]')*mean_Iint' ; % pinv is the pseudo inverse of the matrix because it is not a square matrix
ExcFact = factors(1) ;
InhFact = factors(2) ;

predicted_Iint = ExcFact*mean_Gexc + InhFact*mean_Ginh ; % 

figure
subplot(2,1,1), plot(mean_Gexc)
hold on, plot(mean_Ginh,'r')
title('mean conductances')
legend('exc','inh')
subplot(2,1,2), plot(mean_Iint)
hold on, plot(predicted_Iint,'r')
title('mean Iint')
legend('actual','predicted')

% clear all variables that need to be found seperately for the mean light response section  
clear Idata frequcutoff* Iexc_fft Iinh_fft Iint_fft mean* G* 


% USE EXCFACT AND INHFACT TO PREDICT THE MEAN LIGHT VARIANCE RESPONSE
% separate data by holding potential and keep only section to be analyzed
Idata{1} = EpochCondition(epochCond_num).EpochData.Data(exc_data,firstsectionmean:lastsectionmean) ;  
Idata{2} = EpochCondition(epochCond_num).EpochData.Data(inh_data,firstsectionmean:lastsectionmean) ;
Idata{3} = EpochCondition(epochCond_num).EpochData.Data(int_data,firstsectionmean:lastsectionmean) ;

% adjust exc and inh current data by their offsets  
Idata{1} = Idata{1} - repmat(Idata_offsets{1},1,length(Idata{1})) ;
Idata{2} = Idata{2} - repmat(Idata_offsets{2},1,length(Idata{2})) ;

% low pass filter currents to remove noise from instrumentation
freqcutoff_exc = freqcutoff_exc/(10000/length(Idata{1})) ; % this adjusts the freq cutoff for the length and assumes a sample rate of 10kHz
freqcutoff_inh = freqcutoff_inh/(10000/length(Idata{2})) ;
freqcutoff_int = freqcutoff_int/(10000/length(Idata{3})) ;

Iexc_fft = fft(Idata{1},[],2) ; % take fft of exc signal
Iexc_fft(:,1+freqcutoff_exc:length(Iexc_fft)-freqcutoff_exc) = 0 ; % cut out high frequencies in first and second half of fft
Idata{1} = real(ifft(Iexc_fft,[],2)) ; % inverse fft

Iinh_fft = fft(Idata{2},[],2) ; % take fft of inh signal
Iinh_fft(:,1+freqcutoff_inh:length(Iinh_fft)-freqcutoff_inh) = 0 ; % cut out high frequencies in first and second half of fft
Idata{2} = real(ifft(Iinh_fft,[],2)) ; % inverse fft

Iint_fft = fft(Idata{3},[],2) ; % take fft of inh signal
Iint_fft(:,1+freqcutoff_int:length(Iint_fft)-freqcutoff_int) = 0 ; % cut out high frequencies in first and second half of fft
Idata{3} = real(ifft(Iint_fft,[],2)) ; % inverse fft

% calculate variance of each signal
var_Iexc = mean(var(Idata{1},0,2)) ;
var_Iinh = mean(var(Idata{2},0,2)) ;
var_Iint = mean(var(Idata{3},0,2)) ;

% calculate variance of Gexc and Ginh
var_Gexc = var_Iexc/(exc_drive^2) ; %var_Iexc = var_Gexc * exc_drive^2
var_Ginh = var_Iinh/(inh_drive^2) ;

% Calculate the predicted intermediate variance based on Iexc and Iinh (if the two currents are uncorrelated)
predicted_var_Int = (var_Gexc*(ExcFact^2))+(var_Ginh*(InhFact^2)) ;

predictedVsActual_varIint = predicted_var_Int/var_Iint ; 

% caluculate the power spectrum of the actual currents 
[powerspec_xvalues_Iint, powerspec_Iint] = PowerSpectrumFinder(Idata{3},samplerate) ;
[powerspec_xvalues_Iexc, powerspec_Iexc] = PowerSpectrumFinder(Idata{1},samplerate) ;
[powerspec_xvalues_Iinh, powerspec_Iinh] = PowerSpectrumFinder(Idata{2},samplerate) ;

% calculate the power spectrum of the predicted intermediate current
pred_powerspec_Iint = ((powerspec_Iexc./exc_drive^2)*(ExcFact^2))+((powerspec_Iinh./inh_drive^2)*(InhFact^2)) ;

PredvActual_powerIint = pred_powerspec_Iint./powerspec_Iint ; % ratio of the predicted over the simulated
smoothed_powerspec_ratio = SmoothPowerSpectrum(PredvActual_powerIint(1:length(powerspec_xvalues_Iint)), powerspec_xvalues_Iint, 3, 1) ;% function to smooth powerspec for log scaling


figure
text('units','normalized') %this sets the text in units where (1,1) is upper right
text(.25, 1,['Cell = ' CellInfo_str],'units','normalized')
text(.25, .8, ['inhibitort factor = ' num2str(InhFact,'%10.2e\n')],'units','normalized')
text(.25, .7, ['excitatory factor = ' num2str(ExcFact,'%10.2e\n')],'units','normalized')
text(.25, .5, ['variance Iexc = ' num2str(var_Iexc,'%10.2e\n')],'units','normalized')
text(.25, .4, ['variance Iinh = ' num2str(var_Iinh,'%10.2e\n')],'units','normalized')
text(.25, .3, ['variance Iint = ' num2str(var_Iint,'%10.2e\n')],'units','normalized')
text(.25, .2, ['predicted Resdid Int = ' num2str(predicted_var_Int,'%10.2e\n')],'units','normalized')
text(.25, .1, ['predicted/actual var Int = ' num2str(predictedVsActual_varIint,'%10.2e\n')],'units','normalized')

figure
subplot(1,2,1), plot(powerspec_xvalues_Iint, powerspec_Iint(1:length(powerspec_xvalues_Iint)))
hold on, plot(powerspec_xvalues_Iint, pred_powerspec_Iint(1:length(powerspec_xvalues_Iint)),'r') 
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz) log')
ylabel('power, log')
title('power spectrum of resdual Iint')
legend('actual','predicted')

subplot(1,2,2),plot(powerspec_xvalues_Iint,PredvActual_powerIint(1:length(powerspec_xvalues_Iint)))
hold on, plot(smoothed_powerspec_ratio.Freq,smoothed_powerspec_ratio.PowerSpec,'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz) log')
ylabel('power, log')
title('ratio of predicted/actual power spectrum of resdual Iint')
legend('ratio','smoothed ratio')



