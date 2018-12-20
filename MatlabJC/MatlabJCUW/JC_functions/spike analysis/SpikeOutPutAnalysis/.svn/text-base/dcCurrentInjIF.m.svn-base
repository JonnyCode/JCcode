function ForIgor = dcIF(Input,Parameters,id,A) ;

% get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
for a = 1:length(epochs) ; % for each spike epoch
    [predata{a}, error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data
    [excg{a}, inhg{a}, error] = ITCReadEpochStmGClamp(epochs(a), 0, fp);
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp); % get sampling interval
    SI(a) = SI(a) * 1e-6; % Sampling interval in sec
end

% interpolate all signals so they are pseudo sampled at the same rate
MSI = min(SI) ; % find the signal with the highest sampling rate
SLsec = ((Parameters.PrePnts+Parameters.StmPnts+Parameters.PostPnts)/10000) ; % entire stimulus length in sec
time = [MSI:MSI:SLsec] ; % time vector in sec 
for a = 1:length(epochs) ; % for each epoch
    data(a,:) = interp1([SI(a):SI(a):SLsec],predata{a},time,'linear','extrap') ;  % interpolate the data
    Gexc(a,:) = interp1([SI(a):SI(a):SLsec],excg{a},time,'linear','extrap') ;  % interpolate the data
    Ginh(a,:) = interp1([SI(a):SI(a):SLsec],inhg{a},time,'linear','extrap') ;  % interpolate the data
end

data = data - 10 ; % account for 10mV liquid junction potential

% use integrate and fire model to create spikes
V_trace = LIFmodelG(Gexc(1,:),Ginh(1,:),25,1/MSI) ;

figure
plot(time,data(1,:))
hold on
plot(time,V_trace(1,:),'r')

