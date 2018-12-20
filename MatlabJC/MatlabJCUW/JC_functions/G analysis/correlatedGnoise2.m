function [G_exc,G_inh] = correlatedGnoise2(CellInfo_str,epochCond_num_str,exc_epochs_str, inh_epochs_str, int_epochs_str, firstsection_str,lastsection_str) ;

% This function will analyze conductances recordings and create figures.

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
firstsection = str2num(firstsection_str) ;
if strncmp(lastsection_str,'end',3) ;
    lastsection = length(EpochCondition(epochCond_num).EpochData.Data) ;
else
lastsection = str2num(lastsection_str) ;
end

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

% separate data by holding potential and keep only section to be analyzed
Idata{1} = EpochCondition(epochCond_num).EpochData.Data(exc_data,firstsection:lastsection) ;  
Idata{2} = EpochCondition(epochCond_num).EpochData.Data(inh_data,firstsection:lastsection) ;
Idata{3} = EpochCondition(epochCond_num).EpochData.Data(int_data,firstsection:lastsection) ;

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

% calculate current residual for each trial at each holding potential
residual_Iexc = Idata{1} - repmat(mean_Iexc,numexc,1) ;
residual_Iinh = Idata{2} - repmat(mean_Iinh,numinh,1) ;
residual_Iint = Idata{3} - repmat(mean_Iint,numint,1) ;

% calculate mean of the variance of each residual at each holding potential
var_Iexc = mean(var(residual_Iexc,0,2)) ;
var_Iinh = mean(var(residual_Iinh,0,2)) ;
var_Iint = mean(var(residual_Iint,0,2)) ;


% PLOT FIGURES
% plot exc currents
figure(1), 
subplot(5,1,1), plot([1:length(Idata{1})],Idata{1})
title('Hold at Gaba reversal potential')
legend show

% plot individual inh currents
subplot(5,1,2), plot([1:length(Idata{2})],Idata{2})
title('Hold at Ampa reversal potential')
legend show

% plot intermediate currents
subplot(5,1,3), plot([1:length(Idata{3})],Idata{3})
title('Hold at intermediate potential')
legend show

% plot mean currents
subplot(5,1,4:5), plot(mean_Iexc(1,:))
hold on, plot(mean_Iinh(1,:),'r')
plot(mean_Iint(1,:),'g')
legend('mean Iexc', 'mean Iinh', 'mean Iint')
xlabel('sample points')
ylabel('Current (pA)')

% plot mean currents and residuals
figure(2), 
subplot(2,1,1), plot(mean_Iexc(1,:),'b--')
hold on, plot([1:length(residual_Iexc)],residual_Iexc)
legend('mean')
title('mean and residual Exc currents')
legend show

subplot(2,1,2), plot(mean_Iinh(1,:),'r--')
hold on, plot([1:length(residual_Iinh)],residual_Iinh)
title ('mean and residual Inh currents')
legend('mean')
xlabel('sample points')
ylabel('current (pA)')
legend show

% plot absolute mean amplitude of current against variance of current 
% this may help assess isolation of current
figure(3), 
subplot(4,1,1), plot(mean_Iexc)
title('mean exc current')
subplot(4,1,2), plot(mean_Iinh,'r')
title('mean inh current')
subplot(4,1,3:4), plot(abs(mean_Iexc),var(Idata{1}))
hold on, plot(abs(mean_Iinh),var(Idata{2}),'r')
title ('current amplitude vs current variance')
xlabel('absolute current amplitude (pA)')
ylabel('variance')



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

% calculate the variance of each each conductance
var_Gexc = var_Iexc/(exc_drive^2) ; %var_Iexc = var_Gexc * exc_drive^2
var_Ginh = var_Iinh/(inh_drive^2) ;

% caculate the residual conductances
residual_Gexc = G_exc - repmat(mean_Gexc,numexc,1) ;
residual_Ginh = G_inh - repmat(mean_Ginh,numinh,1) ;

% calculate an intermediate current by finding two factors excFact and
% inhFact such that mean(excFact*Iexc + inhFact*Iinh - Iint)^2 is minimized  

factors = pinv([mean_Gexc;mean_Ginh]')*mean_Iint' ; % pinv is the pseudo inverse of the matrix because it is not a square matrix
ExcFact = factors(1) ;
InhFact = factors(2) ;

predicted_Iint = ExcFact*mean_Gexc + InhFact*mean_Ginh ; % 

predicted2_Iint = (exc_drive/2)*mean_Gexc + (inh_drive/2)*mean_Ginh ; % predicted if you are holding directly inbetween the two holding potentials, you isolated, and they sum linearly 

% calculate the mean squared error of the predicted and actual intermediate

predicted_mse = mean((mean_Iint - predicted_Iint).^2) ; 
predicted2_mse = mean((mean_Iint - predicted2_Iint).^2) ;

mse = mean(((repmat(mean_Iint,numint,1) - Idata{3}).^2),2) ; % mse for actual conductance
mean_mse = mean(mse) ;
 
% Calculate the predicted intermediate variance based Iexc and Iinh if the two currents are uncorrelated

predicted_var_Int = (var_Gexc*(ExcFact^2))+(var_Ginh*(InhFact^2)) ;

predictedVsActual_varIint = predicted_var_Int/var_Iint ; 

% calculate cross correlation of mean exc and inh conductances 
Gcrosscorr = xcov(mean_Gexc, mean_Ginh, 'coeff') ;
xvalues = [1:length(Gcrosscorr)]-((length(Gcrosscorr)+1)/2) ;

% calculate auttocorrelation of mean exc and inh conductances
Gautto_exc = xcorr(mean_Gexc,'coeff') ;
Gautto_inh = xcorr(mean_Ginh,'coeff') ;

% find peak and lag of cross corr
Gcc_peak = max(Gcrosscorr) ;
Gcc_trough = min(Gcrosscorr) ;
Gcc_peaklag = xvalues(find(Gcrosscorr == Gcc_peak)) ; %#ok<FNDSB>
Gcc_troughlag = xvalues(find(Gcrosscorr == Gcc_trough)) ; %#ok<FNDSB>


% create amplitude conductance histograms
hist_Gexc = hist(mean_Gexc,amp_bins) ;  %amplitude histogram exc
hist_Ginh = hist(mean_Ginh,amp_bins) ;

% calculate the mean power spectrum of each conductance
[powerspec_xvalues_Gexc, mean_powerspec_Gexc] = PowerSpectrumFinder(G_exc,samplerate) ;
[powerspec_xvalues_Ginh, mean_powerspec_Ginh] = PowerSpectrumFinder(G_inh,samplerate) ;

% calculate the power spectrum of the mean of exc and inh conductances
[powerspec_xvalues_Gexc, powerspec_mean_Gexc] = PowerSpectrumFinder(mean_Gexc,samplerate) ;
[powerspec_xvalues_Ginh, powerspec_mean_Ginh] = PowerSpectrumFinder(mean_Ginh,samplerate) ;

% calculate the mean power spectrum of the residuals conductances
[powerspec_xvalues_Gexc, powerspec_residual_Gexc] = PowerSpectrumFinder(residual_Gexc,samplerate) ;
[powerspec_xvalues_Ginh, powerspec_residual_Ginh] = PowerSpectrumFinder(residual_Ginh,samplerate) ;

% calculate the difference between...(see below)
powerspec_diff_Gexc = powerspec_mean_Gexc+powerspec_residual_Gexc-mean_powerspec_Gexc ;
powerspec_diff_Ginh = powerspec_mean_Ginh+powerspec_residual_Ginh-mean_powerspec_Ginh ;

% caluculate the power spectrum of the actual residual current int
[powerspec_xvalues_ResIint, powerspec_residual_Iint] = PowerSpectrumFinder(residual_Iint,samplerate) ;

% calculate the predicted power spectrum of the intermediate resiudal
% current assuming indepedance
pred_powerspec_ResIint = (powerspec_residual_Gexc*(ExcFact^2))+(powerspec_residual_Ginh*(InhFact^2)) ;

PredvActual_powerResIint = pred_powerspec_ResIint./powerspec_residual_Iint ; % ratio of the predicted over the simulated
smoothed_powerspec_ratio = SmoothPowerSpectrum(PredvActual_powerResIint(1:length(powerspec_xvalues_ResIint)), powerspec_xvalues_ResIint, 3, 1) ;% function to smooth powerspec for log scaling


% PLOT FIGURES
% plot individual exc conductances
figure(4), 
subplot(2,1,1), plot([1:length(G_exc)],G_exc)
title('exc conductances')
ylabel('conductance (nS)')
legend show

% plot individual inh conductances
subplot(2,1,2), plot([1:length(G_inh)],G_inh)
title('inh conductance')
xlabel('sample points')
ylabel('conductance (nS)')
legend show


% plot the actual intermediate current, best predicted current, and predicted current at hold -30mV 
figure(5), 
plot(predicted_Iint)
hold on
plot(predicted2_Iint,'r')
plot(mean_Iint,'g')
title('intermediate currents')
xlabel('sample points')
ylabel('current (pA)')
legend('best predicted','predicted at -30','actual')


% plot cross correlation of conductances
figure(6), 
subplot(3,1,1), plot(xvalues,Gcrosscorr)
hold on, plot(Gcc_peaklag,Gcc_peak,'k*')
plot(Gcc_troughlag,Gcc_trough, 'k*')
title('cross correlation of conductances')
xlabel('sample points')
ylabel('cc coefficient')

%plot attocorrelations
subplot(3,1,2:3), plot(xvalues, Gautto_exc)
hold on, plot(xvalues, Gautto_inh,'r')
title('auttocorrelation of conductances')
legend('exc', 'inh')

% plot conductances adjusted by lag 
if Gcc_peaklag>0 ;  % if the peak lag is positive ...
    adjusted_mean_Gexc = mean_Gexc(Gcc_peaklag+1:end) ;
    adjusted_mean_Ginh = mean_Ginh(1:end-Gcc_peaklag) ;
else
    adjusted_mean_Gexc = mean_Gexc(1:end+Gcc_peaklag) ;
    adjusted_mean_Ginh = mean_Ginh(-Gcc_peaklag+1:end) ;
end

figure (7)  
subplot(3,1,1),plot(adjusted_mean_Gexc)
hold on, plot(adjusted_mean_Ginh,'r') ;
title('conductances adjusted for lag')
xlabel('sample points')
ylabel('conductance (nS)')

% plot cross correlation of adjusted conductances to assure they were adjusted correctly  
testcc = xcov(adjusted_mean_Gexc,adjusted_mean_Ginh,'coeff') ;
test_xvalues = [1:length(testcc)]-((length(testcc)+1)/2) ;
subplot(3,1,2), plot(test_xvalues, testcc)
title('crosscorrelation of adjusted conductances')
xlabel('sample points')
ylabel('cc coefficient')

% plot the amplitude exc vs amplitude of inh after adjustment for lag
subplot(3,1,3), plot(smooth(adjusted_mean_Gexc,100),smooth(adjusted_mean_Ginh,100),'*') ;
title('adjusted conductances')
xlabel('smoothed amplitude exc')
ylabel('smoothed amplitude inh')


% plot mean amplitude vs. variance of conductances
figure(8), plot(smooth(mean_Gexc,100),smooth(var(G_exc),100),'.')
hold on, plot(smooth(mean_Ginh,100),smooth(var(G_inh),100),'r.')
title ('conductance amplitude vs variance')
xlabel('smoothed conductance amplitude (nS)')
ylabel('smoothed conducance variance')
legend('exc','inh')

% plot amplitude histogram of mean conductance
figure(9), plot(hist_Gexc) ;
hold on, plot(hist_Ginh,'r') ;
title('mean conductance amplitude histogram')
xlabel('conductance amplitude (pS)')
ylabel('number observations')
legend('exc', 'inh')

% plot mean power spectrum of the conductances
figure(10),
subplot(2,2,1), plot(powerspec_xvalues_Gexc,mean_powerspec_Gexc(1:length(powerspec_xvalues_Gexc)))
hold on, plot(powerspec_xvalues_Ginh, mean_powerspec_Ginh(1:length(powerspec_xvalues_Ginh)),'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz) log')
ylabel('power log')
title('mean power spectrum of conductances')
legend('exc', 'inh')

% plot power spectrums of mean conductances
subplot(2,2,2), plot(powerspec_xvalues_Gexc,powerspec_mean_Gexc(1:length(powerspec_xvalues_Gexc)))
hold on, plot(powerspec_xvalues_Ginh, powerspec_mean_Ginh(1:length(powerspec_xvalues_Ginh)),'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz) log')
ylabel('power, log')
title('power spectrum of mean conductance')

% plot mean power spectrums of residual conductances
subplot(2,2,3), plot(powerspec_xvalues_Gexc,powerspec_residual_Gexc(1:length(powerspec_xvalues_Gexc)))
hold on, plot(powerspec_xvalues_Ginh, powerspec_residual_Ginh(1:length(powerspec_xvalues_Ginh)),'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz) log')
ylabel('power, log')
title('mean power spectrum of conductance residuals')

% plot power spectrum of mean + power spectrum of residual - mean power spectrum
subplot(2,2,4), plot(powerspec_xvalues_Gexc,powerspec_diff_Gexc(1:length(powerspec_xvalues_Gexc))) 
hold on, plot(powerspec_xvalues_Ginh,powerspec_diff_Ginh(1:length(powerspec_xvalues_Ginh)),'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz) log')
ylabel('power, log')
title('powerspec mean + residuals - mean powerspec')

% plot power sectrum of actual and predicted intermediate residual current

figure(11)
subplot(1,2,1), plot(powerspec_xvalues_ResIint, powerspec_residual_Iint(1:length(powerspec_xvalues_ResIint)))
hold on, plot(powerspec_xvalues_ResIint, pred_powerspec_ResIint(1:length(powerspec_xvalues_ResIint)),'r') 
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz) log')
ylabel('power, log')
title('power spectrum of resdual Iint')
legend('actual','predicted')

subplot(1,2,2),plot(powerspec_xvalues_ResIint,PredvActual_powerResIint(1:length(powerspec_xvalues_ResIint)))
hold on, plot(smoothed_powerspec_ratio.Freq,smoothed_powerspec_ratio.PowerSpec,'r')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
xlabel('frequency (Hz) log')
ylabel('power, log')
title('ratio of predicted/actual power spectrum of resdual Iint')
legend('ratio','smoothed ratio')

% ACCESS TO CELL
% plot offset of all traces and highlight those analyzed
figure(12)
plot(EpochCondition(epochCond_num).EpochNumbers,EpochCondition(epochCond_num).EpochData.Offset,'*')
hold on, plot([exc_epochs,inh_epochs,int_epochs],EpochCondition(epochCond_num).EpochData.Offset([exc_data,inh_data,int_data]),'r*')
title('Offset currents')
xlabel('epoch number')
ylabel('offset current (pA)')
legend('trials in epochcondition','trials analyzed')


% FOR POSTERITY
% to be printed
figure (13)
subplot(2,2,1:2), plot(mean_Gexc) 
hold on, plot(mean_Ginh, 'r')
title(CellInfo_str)
xlabel('sample points')
ylabel('mean conductance (nS)')
legend('mean exc', 'mean inh')

subplot(2,2,3), plot(xvalues(length(xvalues)/2:(length(xvalues)/2)+20000),Gautto_exc(length(xvalues)/2:(length(xvalues)/2)+20000))
hold on, plot(xvalues(length(xvalues)/2:(length(xvalues)/2)+20000),Gautto_inh(length(xvalues)/2:(length(xvalues)/2)+20000),'r')
title('auttocorrelations (xcorr)')
xlabel('sample points')
ylabel('cc coefficient')
legend('exc','inh')

subplot(2,2,4), plot(xvalues(54000:66000),Gcrosscorr(54000:66000))
% hold on, plot(Gcc_peaklag,Gcc_peak,'k*')
% plot(Gcc_troughlag,Gcc_trough, 'k*')
title('cross correlation of conductances (xcov)')
xlabel('sample points')
ylabel('cc coefficient')

% access resitance indicators and conductance predictions
% plot the actual intermediate current, best predicted current, and predicted current at hold -30mV 
figure(14)
subplot(2,2,3:4), plot(predicted_Iint)
hold on, plot(predicted2_Iint,'r')
hold on, plot(mean_Iint,'g')
title('intermediate currents')
xlabel('sample points')
ylabel('current (pA)')
legend('best predicted','predicted at -30','actual')

% plot offset of all traces and highlight those analyzed
subplot(2,2,2), plot(EpochCondition(epochCond_num).EpochNumbers,EpochCondition(epochCond_num).EpochData.Offset,'*')
hold on, plot([exc_epochs,inh_epochs,int_epochs],EpochCondition(epochCond_num).EpochData.Offset([exc_data,inh_data,int_data]),'r*')
title('Offset currents')
xlabel('epoch number')
ylabel('offset current (pA)')
text('units','normalized') %this sets the text in units where (1,1) is upper right
text(.1, .9,['exc eppochs = ' exc_epochs_str],'units','normalized')
text(.1, .8,['inh eppochs = ' inh_epochs_str],'units','normalized')
text(.1, .7,['int eppochs = ' int_epochs_str],'units','normalized')

% add the following text before printing
subplot(2,2,1), 
text('units','normalized') %this sets the text in units where (1,1) is upper right
text(.25, 1,['Cell = ' CellInfo_str],'units','normalized')
text(.25, .9,['best predicted mse = ' num2str(predicted_mse,'%10.2e\n')],'units','normalized')
text(.25, .8,['predicted mse at -30mV = ' num2str(predicted2_mse,'%10.2e\n')],'units','normalized') 
text(.25, .7,['mean mse = ' num2str(mean_mse,'%10.2e\n')],'units','normalized')

text(.25, .5, ['variance Resid Iexc = ' num2str(var_Iexc,'%10.2e\n')],'units','normalized')
text(.25, .4, ['variance Resid Iinh = ' num2str(var_Iinh,'%10.2e\n')],'units','normalized')
text(.25, .3, ['variance Resid Iint = ' num2str(var_Iint,'%10.2e\n')],'units','normalized')
text(.25, .2, ['predicted var Resdid Int = ' num2str(predicted_var_Int,'%10.2e\n')],'units','normalized')
text(.25, .1, ['predicted/actual var Resid Int = ' num2str(predictedVsActual_varIint,'%10.2e\n')],'units','normalized')

