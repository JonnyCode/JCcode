function ForIgor = BasicGAnalysis(Input,Parameters,id,A) ;

% this function will analysize conductances from rand seed exp, plot the 
% coductances produce the power spectrum and power spectrum of residual traces.

% get data

[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
for a = 1:length(epochs) ; % for each spike epoch
    [predata{a}, error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data
    [prestm{a}, error] = ITCReadEpochStm(epochs(a), 0, fp);  % get the light stimulus
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp); % get sampling interval
    SI(a) = SI(a) * 1e-6; % Sampling interval in sec
end

% interpolate all signals so they are pseudo sampled at the same rate
MSI = min(SI) ; % find the signal with the highest sampling rate
SLsec = ((Parameters.PrePnts+Parameters.StmPnts+Parameters.PostPnts)/10000) ; % entire stimulus length in sec
time = [MSI:MSI:SLsec] ; % time vector in sec 
for a = 1:length(epochs) ; % for each epoch
    data(a,:) = interp1([SI(a):SI(a):SLsec],predata{a},time,'linear','extrap') ;  % interpolate the data
    stm(a,:) = interp1([SI(a):SI(a):SLsec],prestm{a},time,'linear','extrap') ;  % interpolate the data
end
data = data-repmat(mean(data(:,1:Parameters.PrePnts),2),1,length(data)) ; % subtract off prestimulus current

% rectifier and zero stim
negstmPnts = find(stm<0) ;  % find indicies which would be getting a negative stim voltage
stm(negstmPnts) = 0 ;       % make these points zero
stm = stm - repmat(mean(stm(:,Parameters.PrePnts:Parameters.PrePnts+Parameters.StmPnts),2),1,size(stm,2)) ; % subtract off mean of stimulus during time varying stimulus

% currents into coductances
if strcmp(id,'Exc')
    data = data/-61.1 ; % change pA into nS
else
    data = data/61.1 ;
end

% get mean power spectrum of all responses and stim
[powerspec_xvalues, mean_powerspec] = PowerSpectrumFinder(data,1/MSI) ;
[powerspec_xvalues, mean_powerspec] = PowerSpectrumFinder(data,1/MSI) ;


identifier = ['PSXexc',num2str(A)] ;
ForIgor.(identifier) = powerspec_xvaluesEXC ;

identifier = ['PSexc',num2str(A)] ;
ForIgor.(identifier) = mean_powerspecEXC ;

figure
plot(powerspec_xvaluesEXC, mean_powerspecEXC(1:length(powerspec_xvaluesEXC)),'g')
hold on
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;
title(identifier)

identifier = ['PSXinh',num2str(A)] ;
ForIgor.(identifier) = powerspec_xvaluesINH ;

identifier = ['PSinh',num2str(A)] ;
ForIgor.(identifier) = mean_powerspecINH ;

gcf
plot(powerspec_xvaluesINH, mean_powerspecINH(1:length(powerspec_xvaluesEXC)),'r')
