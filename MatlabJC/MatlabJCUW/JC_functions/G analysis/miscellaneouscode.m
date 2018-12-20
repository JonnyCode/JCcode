% miscellaneous set of code previously below Ganalysis (may be usefull but
% I cannot remmeber what I was using it for)


%% Get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

%% Prespike stimulus (STA) from cell attached recording
if perform.STACA == 1 ;

epochs = str2num(Input(A).CA) ;
for a = 1:length(epochs) ; % for each spike epoch
    [stm(a,:), error] = ITCReadEpochStm(epochs(a), 0, fp);  % get the light stimulus
    [data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp);    % get cell attached data
end

data = data - repmat(mean(data(:,1:Parameters.PrePnts),2),1,size(data,2)) ; % subtract mean of prepoints 

% rectifier and zero stim
negstmPnts = find(stm<0) ;  % find indicies which would be getting a negative stim voltage
stm(negstmPnts) = 0 ;       % make these points zero
stm = stm - repmat(mean(stm(:,Parameters.PrePnts:Parameters.PrePnts+Parameters.StmPnts),2),1,size(stm,2)) ; % subtract off mean of stimulus during time varying stimulus

% detect spikes in cell attached data
SpikePnts = SpikeDetection(data,10,10000) ; % data,threshold,samplerate

% get STA
% find spikes that can be used for STA
for a = 1:length(SpikePnts) ;                                                               % for each spike epoch ...
    STA_spikes{a} = find(SpikePnts{a}>Parameters.PrePnts+Parameters.STAPnts & SpikePnts{a}<Parameters.PrePnts+Parameters.StmPnts) ;    % indicies of SpikePnts vector during the stimuli and far enough out to get the prespike wave form 
    NumSpikes(a) = length(STA_spikes{a}) ;                                                  % number of spikes to be used per trial for STA
end % end epoch loop

SumSpikes = [0,cumsum(NumSpikes)] ;

% create a matrix of NaNs for all the pre spike wave forms
PreSpike_stm = nans(SumSpikes(end),floor(Parameters.STAPnts/Parameters.DecimatePnts)) ; % the STA will be decimated, thus the decimate

for a = 1:length(SpikePnts) ;                                                               % for each spike epoch ...
    for b = 1:length(STA_spikes{a}) ; % for each spike that can be used to note a prespike waveform...   
        PreSpike_stm(SumSpikes(a)+b,:) = DecimateWave(stm(a,SpikePnts{a}(STA_spikes{a}(b))-Parameters.STAPnts+1:SpikePnts{a}(STA_spikes{a}(b))),Parameters.DecimatePnts) ;                   
    end % spike loop    
end % end epoch loop

% get average pre spike waveforms over all epochs   
STA = mean(PreSpike_stm) ;
STA = STA/max(abs(STA)) ;    % normalize to maximum 

% prep structure Igor export
identifier = ['STACA',num2str(A)] ;
ForIgor.(identifier) = STA ;

figure
plot(STA)
title(identifier)

clearvars -except A Input Parameters  perform fp ForIgor 
end

%% Linear Filter from light stimuli to Gexc
if perform.LinFExc == 1 ;
    id = 'Exc' ;
    Temp = light2gFilters(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end % end linear filter loop

if perform.LinFInh == 1 ;
    id = 'Inh' ;
    Temp = light2gFilters(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end % end linear filter loop

if perform.LinFInhApb == 1 ;
    id = 'InhApb' ;
    Temp = light2gFilters(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end % end linear filter loop


%% esimate prespike exc and inh synaptic currents from Vclamp and Iclamp data

if perform.STIIC == 1 ;

% get iclamp example cell data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(11).cellname]);

epochsIC = str2num(Input(11).IC) ; % control inh data
for a = 1:length(epochsIC) ; % for each spike epoch
    [stmIC(a,:), error] = ITCReadEpochStm(epochsIC(a), 0, fp);  % get the light stimulus
    [seedIC(a),error]= ITCGetSeed(epochsIC(a),0,0, fp) ; % seed used 
    [dataIC(a,:), error] = ITCReadEpoch(epochsIC(a), 0, fp);    % get cell attached data
    [SIIC(a), error] = ITCGetSamplingInterval(epochsIC(a), fp);
    SRIC(a) = 1/(SIIC(a) * 1e-6); % Sampling rate in Hz
end
dataIC = dataIC - 10 ; % correct for liquid junction potential
stmIC = stmIC - repmat(mean(stmIC,2),1,size(stmIC,2)) ; % subtract off mean of stimulus 
stmIC = stmIC./repmat(max(stmIC,[],2),1,size(stmIC,2)) ; % normalize each stimulus 

% get vclamp cell exc data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

epochsExc = str2num(Input(A).Exc) ; % control inh data
for a = 1:length(epochsExc) ; % for each spike epoch
    [stmExc(a,:), error] = ITCReadEpochStm(epochsExc(a), 0, fp);  % get the light stimulus
    [seedExc(a),error]= ITCGetSeed(epochsExc(a),0,0, fp) ; % seed used 
    [dataExc(a,:), error] = ITCReadEpoch(epochsExc(a), 0, fp);    % get cell attached data
    [SIExc(a), error] = ITCGetSamplingInterval(epochsExc(a), fp);
    SRExc(a) = 1/(SIExc(a) * 1e-6); % Sampling rate in Hz
end
dataExc = dataExc - repmat(mean(dataExc(:,1:Parameters.PrePnts),2),1,size(dataExc,2)) ; % subtract mean of prepoints 
stmExc = stmExc - repmat(mean(stmExc,2),1,size(stmExc,2)) ; % subtract off mean of stimulus 
stmExc = stmExc./repmat(max(stmExc,[],2),1,size(stmExc,2)) ; % normalize each stimulus 

% get vclamp cell inh data
epochsInh = str2num(Input(A).Inh) ; % control inh data
for a = 1:length(epochsInh) ; % for each spike epoch
    [stmInh(a,:), error] = ITCReadEpochStm(epochsInh(a), 0, fp);  % get the light stimulus
    [seedInh(a),error]= ITCGetSeed(epochsInh(a),0,0, fp) ; % seed used 
    [dataInh(a,:), error] = ITCReadEpoch(epochsInh(a), 0, fp);    % get cell attached data
    [SIInh(a), error] = ITCGetSamplingInterval(epochsInh(a), fp);
    SRInh(a) = 1/(SIInh(a) * 1e-6); % Sampling rate in Hz
end
dataInh = dataInh - repmat(mean(dataInh(:,1:Parameters.PrePnts),2),1,size(dataInh,2)) ; % subtract mean of prepoints 
stmInh = stmInh - repmat(mean(stmInh,2),1,size(stmInh,2)) ; % subtract off mean of stimulus 
stmInh = stmInh./repmat(max(stmInh,[],2),1,size(stmInh,2)) ; % normalize each stimulus 

% find unique light seed stimuli used to record current clamped cell
uSet = unique(seedIC) ;

for a = 1:length(uSet) ; % for each unique stimuli used
    IndIC{a} = find(seedIC==uSet(a)) ; % identify the rows where the this stim was used in the iclamp data
    IndExc{a} = find(seedExc==uSet(a)) ; % "" in exc data
    IndInh{a} = find(seedInh==uSet(a)) ; % "" in inh data
    numTrials(a) = min([length(IndIC{a}),length(IndExc{a}),length(IndInh{a})]) ; % the minimum number of repeated trials in exc,inh, or iclamp
end

% assumed reversal potentials
Rexc = 0 ; 
Rinh = -80 ;

% prep matrix
V = nans(sum(numTrials),length(dataIC));
Iexc = V;
Iinh = V;
Gexc = V;
Ginh = V;
round = 0 ;

% estimate exc synaptic currents from synaptic conductances and icalmp data
for a = 1:length(uSet) ; % for each unique stimuli
    for b=1:numTrials(a) ; % for each trial of this unique stim
        round = round+1 ;
        V(round,:) = dataIC(IndIC{a}(b),:) ; % get iclamp data 
        Gexc(round,:) = dataExc(IndExc{a}(b),:)/(Rinh-Rexc) ; % conductances from currents recorded
        Ginh(round,:) = dataInh(IndInh{a}(b),:)/(Rexc-Rinh) ; % conductances from currents recorded
        Iexc(round,:) = Gexc(round,:).*(V(round,:)-Rexc) ; % calculate Iexc and put into matrix
        Iinh(round,:) = Ginh(round,:).*(V(round,:)-Rinh) ;
        Cinh(round,:) = abs(Iinh(round,:))./(abs(Iexc(round,:))+abs(Iinh(round,:))) ; % calculate fraction of inh contribution to total synaptic current
    end
    startpnt = max(cumsum(numTrials(1:a)))-numTrials(a)+1 ;
    endpnt = max(cumsum(numTrials(1:a))) ;
    VMean(a,:)=mean(V(startpnt:endpnt,:)) ;
    GexcMean(a,:)=mean(Gexc(startpnt:endpnt,:)) ;
    GinhMean(a,:)=mean(Ginh(startpnt:endpnt,:)) ;
    IexcMean(a,:)=Gexc(a,:).*(VMean(a,:)-Rexc) ;
    IinhMean(a,:)=Ginh(a,:).*(VMean(a,:)-Rinh) ;
end

% get spike data from I clamp data
SpikePnts = SpikeDetection_WC(V,-20,10000) ; % data,threshold,samplerate

% get STI and STV
% find spikes that can be used for STA
for a = 1:length(SpikePnts) ;                                                               % for each spike epoch ...
    STA_spikes{a} = find(SpikePnts{a}>Parameters.PrePnts+Parameters.STAPnts & SpikePnts{a}<Parameters.PrePnts+Parameters.StmPnts) ;    % indicies of SpikePnts vector during the stimuli and far enough out to get the prespike wave form 
    NumSpikes(a) = length(STA_spikes{a}) ;                                                  % number of spikes to be used per trial for STA
end % end epoch loop

SumSpikes = [0,cumsum(NumSpikes)] ;

% create a matrix of NaNs for all the pre spike wave forms
PreSpike_V = nans(SumSpikes(end),floor(Parameters.STAPnts/Parameters.DecimatePnts)) ; % the STA will be decimated, thus the decimate
PreSpike_Iexc = nans(SumSpikes(end),floor(Parameters.STAPnts/Parameters.DecimatePnts)) ;
PreSpike_Iinh = nans(SumSpikes(end),floor(Parameters.STAPnts/Parameters.DecimatePnts)) ;
PreSpike_Cinh = nans(SumSpikes(end),floor(Parameters.STAPnts/Parameters.DecimatePnts)) ;

for a = 1:length(SpikePnts) ;                                                               % for each spike epoch ...
    for b = 1:length(STA_spikes{a}) ; % for each spike that can be used to note a prespike waveform...   
        PreSpike_V(SumSpikes(a)+b,:) = DecimateWave(V(a,SpikePnts{a}(STA_spikes{a}(b))-Parameters.STAPnts+1:SpikePnts{a}(STA_spikes{a}(b))),Parameters.DecimatePnts) ;   
        PreSpike_Iexc(SumSpikes(a)+b,:) = DecimateWave(Iexc(a,SpikePnts{a}(STA_spikes{a}(b))-Parameters.STAPnts+1:SpikePnts{a}(STA_spikes{a}(b))),Parameters.DecimatePnts) ;
        PreSpike_Iinh(SumSpikes(a)+b,:) = DecimateWave(Iinh(a,SpikePnts{a}(STA_spikes{a}(b))-Parameters.STAPnts+1:SpikePnts{a}(STA_spikes{a}(b))),Parameters.DecimatePnts) ;
        PreSpike_Cinh(SumSpikes(a)+b,:) = DecimateWave(Cinh(a,SpikePnts{a}(STA_spikes{a}(b))-Parameters.STAPnts+1:SpikePnts{a}(STA_spikes{a}(b))),Parameters.DecimatePnts) ;
    end % spike loop    
end % end epoch loop

% Iexc and Iinh distributions time approaches spike time
[hist1exc,bin1exc]=hist(PreSpike_Iexc(:,50),100) ;
[hist2exc,bin2exc]=hist(PreSpike_Iexc(:,200),100) ;
[hist3exc,bin3exc]=hist(PreSpike_Iexc(:,280),100) ;
[hist4exc,bin4exc]=hist(PreSpike_Iexc(:,290),100) ;

[hist1inh,bin1inh]=hist(PreSpike_Iinh(:,50),100) ;
[hist2inh,bin2inh]=hist(PreSpike_Iinh(:,200),100) ;
[hist3inh,bin3inh]=hist(PreSpike_Iinh(:,280),100) ;
[hist4inh,bin4inh]=hist(PreSpike_Iinh(:,290),100) ;


% get average pre spike waveforms over all epochs   
STV = mean(PreSpike_V) ;
STIexc = mean(PreSpike_Iexc) ;
STIinh = mean(PreSpike_Iinh) ;
STCinh = mean(PreSpike_Cinh) ;
STCinh2 = abs(STIinh)./(abs(STIinh)+abs(STIexc)) ;

% prep structure Igor export
identifier = ['STV',num2str(A)] ;
ForIgor.(identifier) = STV ;

figure
plot(STV)
title(identifier)

identifier = ['STIexc',num2str(A)] ;
ForIgor.(identifier) = STIexc ;

figure
plot(STIexc)
title(identifier)

identifier = ['STIinh',num2str(A)] ;
ForIgor.(identifier) = STIinh ;

figure
plot(STIinh)
title(identifier)

identifier = ['STCinh',num2str(A)] ;
ForIgor.(identifier) = STCinh ;

figure
plot(STCinh)
title(identifier)

identifier = ['STCinh2',num2str(A)] ;
ForIgor.(identifier) = STCinh2 ;

figure
plot(STCinh2)
title(identifier)

clearvars -except A Input Parameters  perform fp ForIgor 
end


%%
% power spectrum of exc and inh
if perform.PSg == 1 ;

epochs = str2num(Input(A).Exc) ;
for a = 1:length(epochs) ; % for each spike epoch
    [stm(a,:), error] = ITCReadEpochStm(epochs(a), 0, fp);  % get the light stimulus
    [dataexc(a,:), error] = ITCReadEpoch(epochs(a), 0, fp);    % get cell attached data
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp);
    SR(a) = 1/(SI(a) * 1e-6); % Sampling rate in Hz
end

[powerspec_xvaluesEXC, mean_powerspecEXC] = PowerSpectrumFinder(dataexc,SR(1)) ;
[powerspec_xvaluesLight, mean_powerspecLight] = PowerSpectrumFinder(stm,SR(1)) ;

epochs = str2num(Input(A).Inh) ;
for a = 1:length(epochs) ; % for each spike epoch
    [stm(a,:), error] = ITCReadEpochStm(epochs(a), 0, fp);  % get the light stimulus
    [datainh(a,:), error] = ITCReadEpoch(epochs(a), 0, fp);    % get cell attached data
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp);
    SR(a) = 1/(SI(a) * 1e-6); % Sampling rate in Hz
end

[powerspec_xvaluesINH, mean_powerspecINH] = PowerSpectrumFinder(datainh,SR(1)) ;

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

identifier = ['PSXlight',num2str(A)] ;
ForIgor.(identifier) = powerspec_xvaluesLight ;

identifier = ['PSlight',num2str(A)] ;
ForIgor.(identifier) = mean_powerspecLight ;

gcf
plot(powerspec_xvaluesINH, mean_powerspecINH(1:length(powerspec_xvaluesEXC)),'r')
plot(powerspec_xvaluesLight, mean_powerspecLight(1:length(powerspec_xvaluesLight)),'k')

end
%% Auttocorr inh
if perform.ACInh == 1 ;

epochs = str2num(Input(A).Inh) ;
for a = 1:length(epochs) ; % for each spike epoch
    [stm(a,:), error] = ITCReadEpochStm(epochs(a), 0, fp);  % get the light stimulus
    [data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp);    % get cell attached data
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp);
    SR(a) = 1/(SI(a) * 1e-6); % Sampling rate in Hz
    ac(a,:) = xcov(data(a,:)) ; % auttocorrelation
end

acMean = mean(ac) ;

acXvalues = ([1:length(ac)]-(length(ac)+1)/2)*1000/SR(1) ;

figure
plot(acXvalues,acMean)

end

%% cross correlation of G On inh with Gexc
if perform.CCexcOninh == 1 ; % OFF channel inh isolated via APB

epochsCon = str2num(Input(A).Inh) ; % control inh data
for a = 1:length(epochsCon) ; % for each spike epoch
    [Constm(a,:), error] = ITCReadEpochStm(epochsCon(a), 0, fp);  % get the light stimulus
    [Condata(a,:), error] = ITCReadEpoch(epochsCon(a), 0, fp);    % get cell attached data
    [ConSI(a), error] = ITCGetSamplingInterval(epochsCon(a), fp);
    ConSR(a) = 1/(ConSI(a) * 1e-6); % Sampling rate in Hz
end
Condata = Condata - repmat(mean(Condata(:,1:Parameters.PrePnts),2),1,size(Condata,2)) ; % subtract mean of prepoints 

epochsOff = str2num(Input(A).InhApb) ; % +apb inh data
for a = 1:length(epochsOff) ; % for each spike epoch
    [Offstm(a,:), error] = ITCReadEpochStm(epochsOff(a), 0, fp);  % get the light stimulus
    [Offdata(a,:), error] = ITCReadEpoch(epochsOff(a), 0, fp);    % get cell attached data
    [OffSI(a), error] = ITCGetSamplingInterval(epochsOff(a), fp);
    OffSR(a) = 1/(OffSI(a) * 1e-6); % Sampling rate in Hz
end
Offdata = Offdata - repmat(mean(Offdata(:,1:Parameters.PrePnts),2),1,size(Offdata,2)) ; % subtract mean of prepoints 

overlap = min(size(Offdata,1),size(Condata,1)) ; % number of traces

if ~isequal(Constm(1:overlap,:),Offstm(1:overlap,:)) ; % if all the stim used are in the correct order 
    disp('data out of order?')
end

Inhdata = Condata(1:overlap,:) - Offdata(1:overlap,:) ; % subtract off inh from total inh to find estimates of on inh 
Inhdata = Inhdata/61.6 ; % on inh conductance

epochsExc = str2num(Input(A).Exc) ; % control inh data
for a = 1:length(epochsCon) ; % for each spike epoch
    [Excstm(a,:), error] = ITCReadEpochStm(epochsExc(a), 0, fp);  % get the light stimulus
    [Excdata(a,:), error] = ITCReadEpoch(epochsExc(a), 0, fp);    % get cell attached data
    [ExcSI(a), error] = ITCGetSamplingInterval(epochsExc(a), fp);
    ExcSR(a) = 1/(ExcSI(a) * 1e-6); % Sampling rate in Hz
end
Excdata = Excdata - repmat(mean(Excdata(:,1:Parameters.PrePnts),2),1,size(Excdata,2)) ; % subtract mean of prepoints 
Excdata = Excdata/-61.6 ;

% get cross correlation of G
for a = 1:min(size(Inhdata,1),size(Excdata,1)) ; % for each matching epoch ...
    cc(a,:) = xcov(Excdata(a,Parameters.PrePnts+1:Parameters.PrePnts+Parameters.StmPnts),Inhdata(a,Parameters.PrePnts+1:Parameters.PrePnts+Parameters.StmPnts),'coef') ;
end
ccmean = mean(cc) ;

ccXvalues = [1:length(ccmean)] - (length(ccmean)+1)/2 ;

figure
plot(ccXvalues,ccmean)

end
% end % Cell (A) loop