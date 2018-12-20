
% Get the following outputs (1= get output, 0= don't) 
% Prespike stimulus (STA) from ...
perform_STACA = 0 ;         % cell attached recording
perform_STAIC = 0 ;         % Current clamp
perform_STADC = 0 ;         % dynamic clamp wiht both exc and inh
perform_STADCnoInh = 0 ;    % dynamic clamp with exc but no timevarying inh  
perform_STADCnoExc = 0 ;    % dynamic clamp with inh but no timevarying exc 
perform_STADCnoExcInh = 0 ; % dynamic clamp with neither timevarying exc or inh

% PSTH 
perform_PSTHCA = 0 ;         % cell attached recording
perform_PSTHIC = 0 ;         % Current clamp
perform_PSTHDC = 1 ;         % dynamic clamp wiht both exc and inh
perform_PSTHDCnoInh = 1 ;    % dynamic clamp with exc but no timevarying inh  
perform_PSTHDCnoExc = 1 ;    % dynamic clamp with inh but no timevarying exc 
perform_PSTHDCnoExcInh = 1 ; % dynamic clamp with neither timevarying exc or inh

% Nonlineatiry 
perform_NonlinDC = 1 ;  
perform_NonlinDCnoInh = 1 ;
perform_NonlinDCnoExc = 1 ;

% using the following cells...
for A = [1:2] ;

    
    
%% Cells for analysis

Input(1).cellname = '070908Bc1' ; % good patch, but identity was questioned
Input(1).IC = '[88:102,229:238]' ;
Input(1).DC = '[251:265,305:329]' ; % with synaptic blockers
%Input(1).DC = '[104:143,169:183,214:228]' ; % without blockers
Input(1).DCnoInh = '[266:279]' ;  % with blockers
%Input(1).DCnoInh = '[144:153,184:193]' ;  % without blockers
Input(1).DCnoExc = '[280:294]' ; % with blockers
%Input(1).DCnoExc = '[159:168,194:203]' ;  % without blockers
Input(1).DCnoExcInh = '[295:304]' ;   % with blockers
%Input(1).DCnoExcInh = '[154:158,204:213]' ;   % without blockers
Input(1).DCstm = '[88:92]' ;    % stimulus used when geting G injected in DC

Input(2).cellname = '061908Bc1b' ;  % good patch but dynamic clamp was not optimized yet
Input(2).IC = '[34:43,110:124]' ;
Input(2).DC = '[175:179,195:199,220:224,240:244]' ; % with synaptic blockers
%Input(2).DC = '[50:54,75:79.100:104,135:139,150:154,165:169]' ; % without blockers
Input(2).DCnoInh = '[180:184,200:209,225:229]' ;  % with blockers
%Input(2).DCnoInh = '[55:59,80:84,140:149]' ;  % without blockers
Input(2).DCnoExc = '[185:189,210:214,230:234]' ; % with blockers
%Input(2).DCnoExc = '[60:64,85:89,125:134]' ;  % without blockers
Input(2).DCnoExcInh = '[190:194,235:239]' ;   % with blockers
%Input(2).DCnoExcInh = '[65:74,90:94,155:164]' ;   % without blockers
Input(2).DCstm = '[34:38]' ;    % stimulus used when geting G injected in DC

Input(3).cellname = '061908Bc1' % same cell as above 
Input(3).CA = '[267:281]' ;

% assuming sample rate at 10 kHz
PrePnts = 10000 ;   % number of points before stimuli starts
StmPnts = 60000 ;   % number of points during stimuli
PostPnts = 1000 ;   % number of points after stimuli ends
STA_pnts = 3000 ;   % number of points used for prespike wave forms
DecimatePnts = 10 ; % number of points averaged together in order to downsample prespike waveforms

%% Get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

%% Prespike stimulus (STA) from cell attached recording
if perform_STACA == 1 ;

epochs = str2num(Input(A).CA) ;
for a = 1:length(epochs) ; % for each spike epoch
    [stm(a,:), error] = ITCReadEpochStm(epochs(a), 0, fp);  % get the light stimulus
    [data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp);    % get cell attached data
end

data = data - repmat(mean(data(:,1:PrePnts),2),1,size(data,2)) ; % subtract mean of prepoints 

% rectifier and zero stim
negstmPnts = find(stm<0) ;  % find indicies which would be getting a negative stim voltage
stm(negstmPnts) = 0 ;       % make these points zero
stm = stm - repmat(mean(stm(:,PrePnts:StmPnts),2),1,size(stm,2)) ; % subtract off mean of stimulus during time varying stimulus

% detect spikes in cell attached data
SpikePnts = SpikeDetection(data,10,10000) ; % data,threshold,samplerate

% get STA
% find spikes that can be used for STA
for a = 1:length(SpikePnts) ;                                                               % for each spike epoch ...
    STA_spikes{a} = find(SpikePnts{a}>PrePnts+STA_pnts & SpikePnts{a}<PrePnts+StmPnts) ;    % indicies of SpikePnts vector during the stimuli and far enough out to get the prespike wave form 
    NumSpikes(a) = length(STA_spikes{a}) ;                                                  % number of spikes to be used per trial for STA
end % end epoch loop

SumSpikes = [0,cumsum(NumSpikes)] ;

% create a matrix of NaNs for all the pre spike wave forms
PreSpike_stm = NaNs(SumSpikes(end),floor(STA_pnts/DecimatePnts)) ; % the STA will be decimated, thus the decimate

for a = 1:length(SpikePnts) ;                                                               % for each spike epoch ...
    for b = 1:length(STA_spikes{a}) ; % for each spike that can be used to note a prespike waveform...   
        PreSpike_stm(SumSpikes(a)+b,:) = DecimateWave(stm(a,SpikePnts{a}(STA_spikes{a}(b))-STA_pnts+1:SpikePnts{a}(STA_spikes{a}(b))),DecimatePnts) ;                   
    end % spike loop    
end % end epoch loop

% get average pre spike waveforms over all epochs   
STA = mean(PreSpike_stm) ;
STA = STA/max(abs(STA)) ;    % normalize to maximum 

% prep structure Igor export
identifier = ['STACA',num2str(A)] ;
ForIgor.(identifier) = STA ;

figure(1)
plot(STA)
title('STACA')

clearvars -except ForIgor Input A perform* PrePnts StmPnts PostPnts STA_pnts DecimatePnts fp
end