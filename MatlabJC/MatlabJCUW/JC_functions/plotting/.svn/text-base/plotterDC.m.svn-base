function ploted = plotterDC(Input,Parameters,id,A) ;

% will plot DC V data and then the data with the normalized conductances

% get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
for a = 1:length(epochs) ; % for each spike epoch
    [predata, error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get cell attached data
    [excg, inhg, error] = ITCReadEpochStmGClamp(epochs(a), 0, fp); % get injected G
    [SI, error] = ITCGetSamplingInterval(epochs(a), fp); % get sampling interval
    SI = SI * 1e-6; % Sampling interval in sec
    time = [SI:SI:SI*length(predata)] ; % time vector

figure
subplot(2,1,1)
plot(time,predata)
subplot(2,1,2)
plot(time,excg/max(excg),'b')
hold on
plot(time,inhg/max(inhg),'r')
plot(time,predata/max(abs(predata)),'g')

end

ploted = 1 ;