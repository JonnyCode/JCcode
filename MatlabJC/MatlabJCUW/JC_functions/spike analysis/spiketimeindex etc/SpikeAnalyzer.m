function SpikeStats = SpikeAnalyzer(CellInfo_str, epochCond_num, epochsUwant_str, pre, post, quiet_time) ;

% this function will analyze spikes and create stats and graphs.

% Input = CellInfo file name, epochCondition number, epochs numbers from which you  want to
% detect spikes(vector in quotes), the time during the epoch previous to spike detection and following spike detection, the time period used to define bursts
% output = mostly self explanitory

% JC 9/7/07


% % UNCOMMMENT IF RUNNING AS SCRIPT
% % parameters
% epochsUwant =[10:19] ; %choose the epochs you want
% epochCond_num = 1; %what epochcondtion do you want
% quiet_time = 70;%ms  the time that best defines minimum interburst interval 
% pre=0  ;
% post = 0 ;
% threshold =10 ;

% change strings into numbers
epochsUwant = str2num(epochsUwant_str) ;

% from parameters
cost = 2/quiet_time ; % this will make it less likely that spike distance meteric moves spikes outside of their bursts

% load CellInfo file
cd ~/data_analysis/Index ;
load (CellInfo_str) ;  % parathesies allow it to load that varaible name

% Get appropriate cell structure format
CellInfo = LoadSCIData(CellInfo,1) ;  % for two amps used 

% this function excludes all the epochs but the ones you want
% Input = CellInfo, EpochConditionNumber, vector of the epochs you want
% Output = CellInfo, b = index of epochs in cellinfo
CellInfo = EpochExcluder(CellInfo,epochCond_num,epochsUwant) ; 

% this must come after EpochExcluder
EpochCondition = LoadAndSmoothEpochCondition(CellInfo,5000,1) ;

% plot the epochs you chose
for a =1:length(epochsUwant) ; %for each epoch you want to detect the spike times from
    datamatrix(a,:) = CellInfo.EpochData.Data{epochsUwant(a)+1} ; % put it in a matrix for plotting ease
end
figure(1)
plot([1:length(datamatrix)],datamatrix) ; % plot them all so you can assess best threshold to use


% select threshold
threshold=input('threshold for spike detection?\n') ;  

% plot the detected spikes?
plotdetection = input('plot the detected spikes? (y or n)','s') ;

% This is Gabe's code to detect spike times from Cell Attached data
% Input = CellInfo file name, EpochCondtion(number), pre and post sample
% points to be ingored, spike threshold for detection
% Output = Cell of vector of spike times for each epoch
[SpikeTimeIndex]=GetSpikeTimes_CellAttachedjc(CellInfo,EpochCondition(epochCond_num), pre, post, threshold) ;

% Gabe's code to get PSTH and spiketrains
[SpikeTrain,PSTH] = ConstructPSTH_JC(SpikeTimeIndex,length(datamatrix)) ;

% This function will find spike bursts and analyze them
% Input = SpikeTimeIndex, time inbetween spikes to be considered in a burst 
% Output = for each epoch: number of spikes, number of bursts, average 
% number of spikes per burst, percentage of spikes not in a burst,
% lensdafasf
% average duration of burst, firing rate of first 20ms of burst 
[BurstStats,FirstBurstSpike_time,LastBurstSpike_time] = BurstFinder3(SpikeTimeIndex,quiet_time) ;

% make vectors for plotting the beginning and ends of each burst
BurstStart = zeros(length(FirstBurstSpike_time),length(datamatrix))  ;  % make a zeros vector
BurstEnd = BurstStart ;                                                 
for c = 1:length(FirstBurstSpike_time) ;                                % for each spike train
      BurstStart(c,(FirstBurstSpike_time{c}*10)) = 15 ;                        % put 15 (just smaller than spike label)
      BurstEnd(c,(LastBurstSpike_time{c}*10)) =15 ;
end
  
     
if plotdetection == 'y' ; %if idicated that you wnated to see spikes that were detected
% plot epochs and detected spikes
for b= 1:length(epochsUwant) 
H = figure(b) ;
% set(H,'Position',[-1590,549,2918,509]) ; % for 2 screens
set(H,'Position',[26,660,1239,390]) ; % for 1 screen
plot(CellInfo.EpochData.Data{epochsUwant(b)+1},'k') ;
hold on, plot([zeros(1,pre),SpikeTrain(b,:)]*20,'r') ; % zeros is for graphing only, spike times will still read starting from pre time
plot([zeros(1,pre),BurstStart(b,:)],'b') % plot the spike counted as the first in the burst
plot([zeros(1,pre),BurstEnd(b,:)],'g')   % plot the spike counted as the last in the burst
A = gca ;
set(A,'xlim',([pre, (length(CellInfo.EpochData.Data{epochsUwant(b)+1})-post)])) ; 
end

% quit function if spike detection was not acurrate
contin = input('continue? (y or n)','s') ;
if contin == 'n' ;
    error('quitting')
end

end % plotdection if loop


% plot PSTH smoothed by 50 points (5ms)
figure ;
plot(smooth(PSTH,50)) ;
title('PSTH (5ms moving average)')
xlabel('sample points')
ylabel('av number spikes within 5ms')


% This is the spike distance metric function as assessed by a cumulative
% probability histogram
[DeltaT, tli_spike, tlj_spike, Percent_Pairs_Quantified, DeltaT_Histogram,X_Values,DeltaT_CumProb] = DeltaT_Distribution_From_SCR(SpikeTimeIndex,cost) ;


%plot the spike distance metric cumulative histogram
figure, plot(X_Values,DeltaT_CumProb) ;
H = gca ;
set(H,'Xscale','log') ;
xlabel('Deta T') ;
ylabel('cumulative probability') ;
title('Spike Distance Metric') ;

% find median value of ditribution
[c,d] = min(abs(DeltaT_CumProb-.5)) ; % c=min value, d=index of min value
median_DeltaT = X_Values(d) ; % what detlaT value is that? 

% FOR POSTERITY
SpikeStats = {CellInfo_str, epochCond_num, epochsUwant_str, pre, post, quiet_time,...
    BurstStats(1),BurstStats(2),BurstStats(3),BurstStats(4),BurstStats(5),BurstStats(6)...
    ,BurstStats(7),BurstStats(8),BurstStats(9),BurstStats(10),BurstStats(11),BurstStats(12),...
    BurstStats(13),BurstStats(14),BurstStats(15),BurstStats(16),median_DeltaT} ;

end %end function