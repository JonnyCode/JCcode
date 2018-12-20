function ForIgor = SpkTrigCurr(Input,Parameters,id,A) ;


% set parameters so I don't have to keep writing parameters....
PrePnts = Parameters.PrePnts  ;    % number of points before stimuli starts
StmPnts = Parameters.StmPnts  ;    % number of points during stimuli
PostPnts = Parameters.PostPnts ;    % number of points after stimuli ends
STAPnts = Parameters.STAPnts ;     % number of points used for prespike wave forms
SmoothPnts = Parameters.SmoothPnts ;   % number of points used to smooth spike train to make PSTH

% get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
for a = 1:length(epochs) ; % for each spike epoch
    [stm(a,:), error] = ITCReadEpochStm(epochs(a), 0, fp) ;  % get the injected stimulus
    [data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get cell attached data
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp); % get sampling interval
    SI(a) = SI(a) * 1e-6; % Sampling interval in sec
end

% zero current stim
stm = stm-repmat(mean(stm(:,PrePnts+STAPnts:PrePnts+StmPnts),2),1,length(stm)) ;

% detect spikes in whole cell attached data
SpikePnts = SpikeDetection_WC(data,-20,10000) ; % data,threshold,samplerate

figure
plot(data(2,:))
hold on
plot(SpikePnts{2},-20,'*r')

% get PSTH and SpikeTrain
[PSTH,SpikeTrain] = PSTHmaker(SpikePnts,SmoothPnts,length(data),1/SI(1)) ;   % PSTHmaker(SpikePnts,SmoothPnts,RawDataLength,SampleRate)

% get prespike current (STI)
STI = ifft(mean(fft(SpikeTrain(:,PrePnts+STAPnts:PrePnts+StmPnts),[],2).*conj(fft(stm(:,PrePnts+STAPnts:PrePnts+StmPnts),[],2)))) ; % sta is the correlation of spike train and stim

figure
plot(-[SI(1):SI(1):length(STI)*SI(1)],STI)
title('STI')

% prep ForIgor struct

identifier = ['STI',id,num2str(A)] ;
ForIgor.(identifier) = STI ;

identifier = ['STInorm',id,num2str(A)] ;
ForIgor.(identifier) = STI/max(STI) ;


