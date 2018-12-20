function ForIgor = NNanalysis(Input,Parameters,id,A) ;

% this script will do a nearest neighbor analysis to assess burstiness 

% get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
for a = 1:length(epochs) ; % for each spike epoch
    [predata{a}, error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get cell attached data
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp); % get sampling interval
    SI(a) = SI(a) * 1e-6; % Sampling interval in sec
end

% interpolate all signals so they are pseudo sampled at the same rate
MSI = min(SI) ; % find the signal with the highest sampling rate
SLsec = ((Parameters.PrePnts+Parameters.StmPnts+Parameters.PostPnts)/10000) ; % entire stimulus length in sec
time = [MSI:MSI:SLsec] ; % time vector in sec 
for a = 1:length(epochs) ; % for each epoch
    data(a,:) = interp1([SI(a):SI(a):SLsec],predata{a},time,'linear','extrap') ;  % interpolate the data
end
if ~strcmp(id,'CA') ; % if it is a whole cell recording
    data = data - 10 ; % account for 10mV liquid junction potential
end

% detect spikes in data
if strcmp(id,'CA') ;
    SpikePnts = SpikeDetection(data,15,(1/MSI)) ; % data,threshold,samplerate
else
    SpikePnts = SpikeDetection_WC(data,Parameters.WCSpikeThresh,(1/MSI)) ; % data,threshold,samplerate
end
    
% change points to be appropriate for new sample rate
PrePnts = round(Parameters.PrePnts*.0001/MSI) ;
StmPnts = round(Parameters.StmPnts*.0001/MSI) ;
PostPnts = round(Parameters.PostPnts*.0001/MSI) ;
STAPnts = round(Parameters.STAPnts*.0001/MSI) ;
DecimatePnts = round(Parameters.DecimatePnts*.0001/MSI) ;
QuietPnts = round(Parameters.QuietPnts*.0001/MSI) ; 


% find spikes during stim only and get total isi
isi = [];
for a = 1:length(SpikePnts) ;                                                               % for each spike epoch ...
    Stm_spikes{a} = find(SpikePnts{a}>PrePnts+STAPnts & SpikePnts{a}<PrePnts+StmPnts) ;     % indicies of SpikePnts vector during the stimuli and far enough out to get the prespike wave form 
    NumStmSpikes(a) = length(Stm_spikes{a}) ;                                                   % number of spikes to be used per trial for STA
    isi=[isi,(MSI*1000)*diff(SpikePnts{a}(Stm_spikes{a}))] ;                                             
end % end epoch loop

% change isi into nearest neighbor
for a= 1:length(isi)-1 ; % for each isi but the last
    tempNN(a) = min(isi(a),isi(a+1)) ;% pick the shorter isi
end
    
NN = [isi(1),tempNN,isi(end)] ; % the 1st and last spikes have no choice

% make cumulative histograms
[NNdis,NNbins] = hist(NN,[0:.5:1000]) ;
NNcumdis = cumsum(NNdis) ;
NNcumdisNorm = NNcumdis/max(NNcumdis) ;

NNcumdisRev = cumsum(fliplr(NNdis)) ;
NNcumdisRevNorm = NNcumdisRev/max(NNcumdisRev) ;
NNbinsRev = fliplr(NNbins) ;

% prep matricies for export
identifier = ['NNcd',id,num2str(A)] ;
ForIgor.(identifier) = NNcumdis ;

identifier = ['NNbins',id,num2str(A)] ;
ForIgor.(identifier) = NNbins ;

identifier = ['NNcdN',id,num2str(A)] ;
ForIgor.(identifier) = NNcumdisNorm ;

identifier = ['NNcdR',id,num2str(A)] ;
ForIgor.(identifier) = NNcumdisRev ;

identifier = ['NNbinsR',id,num2str(A)] ;
ForIgor.(identifier) = NNbinsRev ;

identifier = ['NNcdRN',id,num2str(A)] ;
ForIgor.(identifier) = NNcumdisRevNorm ;


% figures
figure
plot(NNbins,NNcumdis)
title(['id',num2str(A)])

figure
plot(NNbins,NNcumdisNorm)
title(['id',num2str(A)])

figure
plot(NNbinsRev,NNcumdisRev)
title(['id',num2str(A)])

figure
plot(NNbinsRev,NNcumdisRevNorm)
title(['id',num2str(A)])

end













    
    
    


