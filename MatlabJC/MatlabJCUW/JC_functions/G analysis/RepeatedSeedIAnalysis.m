function ForIgor = RepeatedSeedIAnalysis(Input,Parameters,id,A) ;

% this function will look at current responses to a repeated light stimulus
% and assess the characterstics of the noise and signal

% 5/4/10 JC

% get data
[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
for a = 1:length(epochs) ; % for each spike epoch
    [stm(a,:), error] = ITCReadEpochStm(epochs(a), 0, fp);  % get the light stimulus
    [data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp);    % get current data
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp); % get sampling interval
    SI(a) = SI(a) * 1e-6; % Sampling interval in sec
end

if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end

[prePnts, error] = ITCGetStmPrePts(epochs(1), 0, 0, fp) ; % points collected beyond which data is not worth analyzing
[postPnts, error] = ITCGetStmTailPts(epochs(1), 0, 0, fp) ;

stm_data = data(:,prePnts:end-postPnts) ; % get only data during stm

% subtract mean off each stm data trial
for a=1:size(stm_data,1) ;
    stm_data(a,:) = stm_data(a,:) - mean(stm_data(a,:)) ;
end

stm_data_mean = mean(stm_data,1) ; % time dependant mean
stm_data_var = var(stm_data,[],1) ; % time dependant variance

ac_mean = xcov(stm_data_mean,'coef') ;
cc_meanvar = xcov(stm_data_mean,stm_data_var,'coef') ;
cc_x = SI(1)*([1:length(cc_meanvar)] - (length(cc_meanvar)+1)/2) ;

ac_mean_halfwidth = abs(cc_x(find(ac_mean>.5,1)))/SI(1) ; % pnts


% select events
for a=1:size(stm_data,1) ;
    eventGenerator(a,:) = smooth(stm_data(a,:).*stm_data_mean,ac_mean_halfwidth) ;

end


time = [1:length(data1)]*SI(1) ;


