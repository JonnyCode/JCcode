% PARAMETERS
load 043008Bc2_a
epochsUwant = [85:164] ;
pre = 0 ;
post = 0 ;
threshold = 10 ;
epochCond_num = 1;
colorA = 'r+' ;
colorB = 'b*' ;
colorC = 'b' ;



% Get appropriate cell structure format
CellInfo = LoadSCIData(CellInfo,1) ;  % for two amps used 

% this function excludes all the epochs but the ones you want
% Input = CellInfo, EpochConditionNumber, vector of the epochs you want
% Output = CellInfo, b = index of epochs in cellinfo
CellInfo = EpochExcluder(CellInfo,epochCond_num,epochsUwant) ; 

% this must come after EpochExcluder
EpochCondition = LoadAndSmoothEpochCondition(CellInfo,5000,1) ;

% This is Gabe's code to detect spike times from Cell Attached data
% Input = CellInfo file name, EpochCondtion(number), pre and post sample
% points to be ingored, spike threshold for detection
% Output = Cell of vector of spike times for each epoch
[SpikeTimeIndex]=GetSpikeTimes_CellAttachedjc(CellInfo,EpochCondition(epochCond_num), pre, post, threshold) ;

for a = 1:8 % each set from a common step contrast

% Gabe's code to get PSTH and spiketrains
[SpikeTrain{a},PSTH{a}] = ConstructPSTH_JC((SpikeTimeIndex(1,a:8:length(epochsUwant))),20000) ;

PSTH_threshfactor = 5 ;

PSTH_smooth{a} = smooth(PSTH{a},100) ;

PSTHsmooth_peak(a) = max(PSTH_smooth{a}(5000:10000)) ;

PSTHsmooth_std(a) = std(PSTH_smooth{a}(1:5000)) ;
PSTH_thresh(a) = PSTH_threshfactor*PSTHsmooth_std(a) ;
PSTHsmooth_threshPnt_start(a) = find(PSTH_smooth{a}>PSTH_thresh(a), 1,'first') ;
PSTHsmooth_threshPnt_end(a) = find(PSTH_smooth{a}>PSTH_thresh(a), 1,'last') ;

spikenumberChange(:,a) = sum(SpikeTrain{a}(:,5000:10000),2) - sum(SpikeTrain{a}(:,1:5000),2) ;


for b = 1:size(SpikeTrain{a},1)
    IndividPSTH{a}(b,:) = smooth(SpikeTrain{a}(b,:),500) ;                                                      % The individuals are smoothed by 50ms
    [IndividPSTH_peak(b,a),IndividPSTH_peakPnt(b,a)] = max(IndividPSTH{a}(b,500:end)) ;                         % find the max of the PSTH and time of max
    IndividPSTH_peakPnt(b,a) = IndividPSTH_peakPnt(b,a) + 500 - 1;
    IndividPSTH_NoiseStd(b,a) = std(IndividPSTH{a}(b,500:5000)) ;                                               % std of prestim data
    IndividPSTH_thresh(b,a) = (IndividPSTH_NoiseStd(b,a) * PSTH_threshfactor) + .0001 ;                                   % threshold = amount above noise cell would believe in signal
    if IndividPSTH_peak(b,a)>IndividPSTH_thresh(b,a) ;                                                          % if there is a point on the PSTH above threshold value ...
    IndividPSTH_threshstart(b,a) = find(IndividPSTH{a}(b,500:end)>IndividPSTH_thresh(b,a), 1,'first') + 500 ;   % find the first point above threshold
    IndividPSTH_threshend(b,a) = find(IndividPSTH{a}(b,IndividPSTH_peakPnt(b,a):end)<IndividPSTH_thresh(b,a), 1,'first') + IndividPSTH_peakPnt(b,a) ;      % find the first point that drops below threshold after the peak 
    else
    IndividPSTH_peak(b,a) = NaN ;
    IndividPSTH_peakPnt = NaN ;
    IndividPSTH_threshstart(b,a) = NaN ;                                                                        % if you cannot find the peak above the noise than you fail to identify light step (put a NaN)
    IndividPSTH_threshend(b,a) = NaN ;    
    end
    
    if sum(SpikeTrain{a}(b,5000:10000))>0 ;
    FirstSpikePnt(b,a) = find(SpikeTrain{a}(b,5000:10000)==1, 1, 'first') + 5000 ;
    LastSpikePnt(b,a) = find(SpikeTrain{a}(b,5000:10000)==1, 1, 'last') + 5000 ;    
    else
    FirstSpikePnt(b,a) = NaN ;    
    LastSpikePnt(b,a) = NaN ;
    end
end

IndividPSTH_mean{a} = mean(IndividPSTH{a}) ;

figure(1)
subplot(5,2,a), plot(IndividPSTH_mean{a},colorC)
hold on
subplot(5,2,9:10), plot(IndividPSTH_mean{a},colorC)
hold on




% subplot(5,2,a), plot(PSTH_smooth{a})
% hold on
% plot(PSTHsmooth_threshPnt_start(a),PSTH_thresh(a),'r*')
% plot(PSTHsmooth_threshPnt_end(a),PSTH_thresh(a),'g*')
% 
% subplot(5,2,9:10), plot(PSTH_smooth{a})
% hold on
% plot(PSTHsmooth_threshPnt_start(a),PSTH_thresh(a),'r*')
% plot(PSTHsmooth_threshPnt_end(a),PSTH_thresh(a),'g*')
end

IndividPSTHpeak_mean = nanmean(IndividPSTH_peak) ;

IndividPSTHthreshstart_mean = nanmean(IndividPSTH_threshstart) ;
IndividPSTHthreshend_mean = nanmean(IndividPSTH_threshend) ;

IndividPSTHduration = IndividPSTH_threshend - IndividPSTH_threshstart ;
IndividPSTHduration_mean = nanmean(IndividPSTHduration) ;

Failures = sum(isnan(IndividPSTH_peak),1) ; % for every trace where the PSTH peak was below the noise threshold

% FirstSpikePnt_duration = LastSpikePnt-FirstSpikePnt ;
% FIrstSpikePntduration_mean = nanmean(FirstSpikePnt_duration) ;
% 
% PSTHsmooth_burstduration = PSTHsmooth_threshPnt_end - PSTHsmooth_threshPnt_start ;

spikenumberChange_mean = nanmean(spikenumberChange) ;

FirstSpikePnt_mean = nanmean(FirstSpikePnt) ;


figure(2)
subplot(3,2,1)
hold on
plot([.125,-.125,.25,-.25,.5,-.5,1,-1], IndividPSTH_peak,colorB)
plot([.125,-.125,.25,-.25,.5,-.5,1,-1],IndividPSTHpeak_mean,colorA)
title('IndividPSTH peak')

subplot(3,2,2)
hold on
plot([.125,-.125,.25,-.25,.5,-.5,1,-1], spikenumberChange,colorB)
plot([.125,-.125,.25,-.25,.5,-.5,1,-1],spikenumberChange_mean,colorA)
title('spikenumberChange')

subplot(3,2,3)
hold on
plot([.125,-.125,.25,-.25,.5,-.5,1,-1], IndividPSTHduration,colorB)
plot([.125,-.125,.25,-.25,.5,-.5,1,-1],IndividPSTHduration_mean,colorA)
title('IndividPSTHduration')

subplot(3,2,4)
hold on
plot([.125,-.125,.25,-.25,.5,-.5,1,-1], IndividPSTH_threshstart,colorB)
plot([.125,-.125,.25,-.25,.5,-.5,1,-1],IndividPSTHthreshstart_mean,colorA)
title('IndividPSTH threshstart')

subplot(3,2,5)
hold on
plot([.125,-.125,.25,-.25,.5,-.5,1,-1], FirstSpikePnt,colorB)
plot([.125,-.125,.25,-.25,.5,-.5,1,-1],FirstSpikePnt_mean,colorA)
title('FirstSpikePnt')

subplot(3,2,6)
hold on
plot([.125,-.125,.25,-.25,.5,-.5,1,-1], Failures,colorB)
title('Failed to reach threshold')


% plot a few detected spikes ontop of real spikes to check detection worked
[SpikeTrainCheck,PSTHCheck] = ConstructPSTH_JC((SpikeTimeIndex),20000) ;
for b= 1:7:length(epochsUwant) 
% set(H,'Position',[-1590,549,2918,509]) ; % for 2 screens
figure
plot(CellInfo.EpochData.Data{epochsUwant(b)+1},'k') ;
hold on, plot([zeros(1,pre),SpikeTrainCheck(b,:)]*20,'r') ; % zeros is for graphing only, spike times will still read starting from pre time
A = gca ;
set(A,'xlim',([pre, (length(CellInfo.EpochData.Data{epochsUwant(b)+1})-post)])) ; 
end

% figure
% plot([.125,-.125,.25,-.25,.5,-.5,1,-1],PSTHsmooth_peak,'k*')
% 
% figure
% plot([.125,-.125,.25,-.25,.5,-.5,1,-1],PSTHsmooth_threshPnt,'k*')
% 
% figure
% plot([.125,-.125,.25,-.25,.5,-.5,1,-1],PSTHsmooth_burstduration,'k*')
% 
% figure
% plot([.125,-.125,.25,-.25,.5,-.5,1,-1],FIrstSpikePntduration_mean,'r+')
% hold on
% plot([.125,-.125,.25,-.25,.5,-.5,1,-1], FirstSpikePnt_duration,'b.')







