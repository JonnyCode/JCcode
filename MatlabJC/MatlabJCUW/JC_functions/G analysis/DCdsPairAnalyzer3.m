function ForIgor = DCdsPairAnalyzer3(Input,Parameters,id,A) ;

% modified from DCdsPairAnalzer2.m took correlation from peak of cross
% correlation (as done in DCdsPairAnalyzer.m not DCdsPairAnalyzer2.m) and added distributed matlab coding 

% JC 6/28/11 
% edit 7/8/11 reverted to making linear analytic predictions with covariance of conductances averaged over
% time (as done in DCdsPairAnalyzer2.m not DCdsPairAnalyzer.m)

% JC 8/13/11 editted to use both cc peak and conductance average residuals.

%% impact of changing pairwise e,i correlations
InhRev = -60 ; % reversal potential in igor code
residualOption = Parameters.residualOption ; % residual option
minusPairwiseShift = 5 ; % number of trials to shift to get rid of pairwise correlations
prePnts = 0 ; % points to avoid before spike data
numBars = 3 ; % number of bars presented in this set of conductances
NullSpikeTime = .2 ; %sec time at the begining of an epoch that spikes will be deleted and ignored in analysis

epochs = str2num(Input(A).(id)) ; %sim and shuffled g

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);

% get data
[epochSize, error] = ITCGetEpochSize(epochs(1), fp) ; % preallocate variables to improve memory usage and speed 
data = nan(length(epochs),epochSize) ;
excg = nan(length(epochs),epochSize) ;
inhg = nan(length(epochs),epochSize) ;

for a = 1:length(epochs) ;
    if Input(A).gHWflag==1 ; % if this is using dynamic clamp through ITC hardware
        [data(a,:), error] = ITCReadEpoch(epochs(a), 2, fp) ; 
        [excg(a,:), inhg(a,:), error] = ITCReadEpochStmGClamp(epochs(a), 3, fp) ;    
    else
        [data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ; 
        [excg(a,:), inhg(a,:), error] = ITCReadEpochStmGClamp(epochs(a), 0, fp) ;
    end
end

% get sampling interval 
[SI, error] = ITCGetSamplingInterval(epochs(1), fp);    
SI = SI * 1e-6; % Sampling interval in sec
if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end
time = [SI:SI:SI*length(data)] ;

% voltage mean
Voltage_mean = mean(data(:)) ;

% spike detection
SpikePnts = SpikeDetection_WC(data,-30,1/SI(1)) ; 

% get rid of spikes within NullSpikeTime
NullSpikePnts = NullSpikeTime/SI(1) ;
for a = 1:length(epochs) ;
    SpikePnts{a}(SpikePnts{a}<=NullSpikePnts) = [] ;
end

numTrials = length(SpikePnts)/4 ; % number of sim or shuffled mimiced trials per cell

% spike thresholds
SpikeThreshold_search = 2 ; 
SpikeThreshold = SpikeThresholdFinder(data,SpikePnts,SI(1),SpikeThreshold_search) ;

SpikeThreshold1_sim = [] ;
SpikeThreshold2_sim = [] ;
SpikeThreshold1_shuff = [] ;
SpikeThreshold2_shuff = [] ;
if SpikeThreshold_search == 1 ;
    for a=1:numTrials ;
        SpikeThreshold1_sim = [SpikeThreshold1_sim,SpikeThreshold{a*4-3}] ;
        SpikeThreshold2_sim = [SpikeThreshold2_sim,SpikeThreshold{a*4-2}] ;
        SpikeThreshold1_shuff = [SpikeThreshold1_shuff,SpikeThreshold{a*4-1}] ;
        SpikeThreshold2_shuff = [SpikeThreshold2_shuff,SpikeThreshold{a*4}] ;
    end
elseif SpikeThreshold_search == 2 ;
    for a=1:numTrials ;
        SpikeThreshold1_sim = [SpikeThreshold1_sim,SpikeThreshold(a*4-3)] ;
        SpikeThreshold2_sim = [SpikeThreshold2_sim,SpikeThreshold(a*4-2)] ;
        SpikeThreshold1_shuff = [SpikeThreshold1_shuff,SpikeThreshold(a*4-1)] ;
        SpikeThreshold2_shuff = [SpikeThreshold2_shuff,SpikeThreshold(a*4)] ;
    end
end
SpikeThreshold_range = [-80:0] ;
SpikeThreshold1_sim_Hist = hist(SpikeThreshold1_sim,SpikeThreshold_range) ;
SpikeThreshold2_sim_Hist = hist(SpikeThreshold2_sim,SpikeThreshold_range) ;
SpikeThreshold1_shuff_Hist = hist(SpikeThreshold1_shuff,SpikeThreshold_range) ;
SpikeThreshold2_shuff_Hist = hist(SpikeThreshold2_shuff,SpikeThreshold_range) ;

SpikeThreshold_Hist = sum([SpikeThreshold1_sim_Hist;SpikeThreshold2_sim_Hist;SpikeThreshold1_shuff_Hist;SpikeThreshold2_shuff_Hist]) ;
SpikeThreshold_mean = sum(SpikeThreshold_Hist.*SpikeThreshold_range)/sum(SpikeThreshold_Hist) ;

% absolute referactory period
absRefPnts = length(data) ; % set an initial absolute referactory period very long
for a=1:length(SpikePnts) ; % for every spike epoch
    if length(SpikePnts{a})>1 ; % if there are more than one spike
       absRefPnts = min([absRefPnts, diff(SpikePnts{a})]) ; % fewest points inbetween neighboring spikes
    end
end
absRef = absRefPnts*SI(1) ; % seconds

% spike trains
spikeTrain1_sim = zeros(numTrials,length(data)) ; 
spikeTrain2_sim = zeros(numTrials,length(data)) ; 
spikeTrain1_shuff = zeros(numTrials,length(data)) ;  
spikeTrain2_shuff = zeros(numTrials,length(data)) ; 

for a = 1:numTrials ;
    spikeTrain1_sim(a,SpikePnts{a*4-3})=1 ; % mimic cell 1 + converging
    spikeTrain2_sim(a,SpikePnts{a*4-2})=1 ; % mimic cell 2 + converging 
    spikeTrain1_shuff(a,SpikePnts{a*4-1})=1 ; % mimic cell 1 - converging
    spikeTrain2_shuff(a,SpikePnts{a*4})=1 ; % mimic cell 2 - converging
end

% spike trains divided by bar
barPnts = floor(length(data)/numBars) ; 

for a=1:numBars ; % preallocate variables to improve memory usage and speed
    spikeTrain1_sim_bar{a} = nan(numTrials,barPnts) ;
    spikeTrain2_sim_bar{a} = nan(numTrials,barPnts) ;
    spikeTrain1_shuff_bar{a} = nan(numTrials,barPnts) ;
    spikeTrain2_shuff_bar{a} = nan(numTrials,barPnts) ;
end

for a=1:numBars ;
    spikeTrain1_sim_bar{a} = spikeTrain1_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    spikeTrain2_sim_bar{a} = spikeTrain2_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    spikeTrain1_shuff_bar{a} = spikeTrain1_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    spikeTrain2_shuff_bar{a} = spikeTrain2_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
end

% spike number for each bar
sn1_sim = nan(numTrials,numBars) ; % preallocate variables to improve memory usage and speed
sn2_sim = nan(numTrials,numBars) ;
sn1_shuff = nan(numTrials,numBars) ;
sn2_shuff = nan(numTrials,numBars) ;
sn1_sim_mean = nan(1,numBars) ;
sn2_sim_mean = nan(1,numBars) ;
sn1_shuff_mean = nan(1,numBars) ;
sn2_shuff_mean = nan(1,numBars) ;

for a=1:numBars ;
    sn1_sim(:,a) = sum(spikeTrain1_sim_bar{a},2) ;
    sn2_sim(:,a) = sum(spikeTrain2_sim_bar{a},2) ;
    sn1_shuff(:,a) = sum(spikeTrain1_shuff_bar{a},2) ;
    sn2_shuff(:,a) = sum(spikeTrain2_shuff_bar{a},2) ;

    sn1_sim_mean(a) = mean(sn1_sim(:,a)) ;
    sn2_sim_mean(a) = mean(sn2_sim(:,a)) ;
    sn1_shuff_mean(a) = mean(sn1_shuff(:,a)) ;
    sn2_shuff_mean(a) = mean(sn2_shuff(:,a)) ;
end

% spike number residuals 
sn1_sim_res = nan(numTrials-2,numBars) ; % preallocate variables to improve memory usage and speed
sn2_sim_res = nan(numTrials-2,numBars) ;
sn1_shuff_res = nan(numTrials-2,numBars) ;
sn2_shuff_res = nan(numTrials-2,numBars) ;

if residualOption == 0 ; % if full mean residual
     for a = 1:numTrials ;
        sn1_sim_res(a,:) = sn1_sim(a,:) - mean(sn1_sim) ;
        sn2_sim_res(a,:) = sn2_sim(a,:) - mean(sn2_sim);
        sn1_shuff_res(a,:) = sn1_shuff(a,:) - mean(sn1_shuff) ;
        sn2_shuff_res(a,:) = sn2_shuff(a,:) - mean(sn2_shuff) ;
     end

elseif residualOption == 1 ; % if local residual calc
    for a = 2:numTrials-1 ;
        sn1_sim_res(a-1,:) = sn1_sim(a,:) - (sn1_sim(a-1,:)+sn1_sim(a+1,:))/2 ;
        sn2_sim_res(a-1,:) = sn2_sim(a,:) - (sn2_sim(a-1,:)+sn2_sim(a+1,:))/2 ;
        sn1_shuff_res(a-1,:) = sn1_shuff(a,:) - (sn1_shuff(a-1,:)+sn1_shuff(a+1,:))/2 ;
        sn2_shuff_res(a-1,:) = sn2_shuff(a,:) - (sn2_shuff(a-1,:)+sn2_shuff(a+1,:))/2 ;
    end
end


% spike number correlation (linear) for each bar
lincoefsigMin = .05 ;
for a=1:numBars ;
    
    % + converging + pairwise 
    [Corr,p,lb,ub] =  corrcoef(sn1_sim_res(:,a),sn2_sim_res(:,a)) ; % correlations 
    sn_plusConverg_plusPairwise_lincoef(a) = Corr(1,2) ;
%    sn_plusConverg_plusPairwise_lincoefsig(a) = p(1,2) ;
%     if sn_plusConverg_plusPairwise_lincoefsig(a)>lincoefsigMin ;
%         sn_plusConverg_plusPairwise_lincoef(a) = 0 ;
%     end
%    sn_plusConverg_plusPairwise_std(a) = (Corr(1,2)-lb(1,2))*sqrt(length(sn1_sim_res(:,a)))/1.96 ; % std from confidence intervals
    Covariation =  cov(sn1_sim_res(:,a),sn2_sim_res(:,a)) ; % covariation 
    sn_plusConverg_plusPairwise_cov(a) = Covariation(1,2) ;
    h=error_ellipse(Covariation,[sn1_sim_mean(a),sn2_sim_mean(a)],'conf',.95) ; % confidence elipse
    sn_plusConverg_plusPairwise_ellipseX(a,:) = get(h,'XData') ;
    sn_plusConverg_plusPairwise_ellipseY(a,:) = get(h,'YData') ;
    sn1_sim_var(a) = Covariation(1,1) ; % variance
    sn2_sim_var(a) = Covariation(2,2) ;
    
    % - converging + pairwise
    [Corr,p,lb,ub] =  corrcoef(sn1_shuff_res(:,a),sn2_shuff_res(:,a)) ; % correlations
    sn_minusConverg_plusPairwise_lincoef(a) = Corr(1,2) ;
%    sn_minusConverg_plusPairwise_lincoefsig(a) = p(1,2) ;
%     if sn_minusConverg_plusPairwise_lincoefsig(a)>lincoefsigMin ;
%         sn_minusConverg_plusPairwise_lincoef(a) = 0 ;
%     end
%    sn_minusConverg_plusPairwise_std(a) = (Corr(1,2)-lb(1,2))*sqrt(length(sn1_sim_res(:,a)))/1.96 ; % std from confidence intervals
    Covariation =  cov(sn1_shuff_res(:,a),sn2_shuff_res(:,a)) ; % covariation
    sn_minusConverg_plusPairwise_cov(a) = Covariation(1,2) ;
    h=error_ellipse(Covariation,[sn1_shuff_mean(a),sn2_shuff_mean(a)],'conf',.95) ; % confidence elipse
    sn_minusConverg_plusPairwise_ellipseX(a,:) = get(h,'XData') ;
    sn_minusConverg_plusPairwise_ellipseY(a,:) = get(h,'YData') ;
    sn1_shuff_var(a) = Covariation(1,1) ; % variance
    sn2_shuff_var(a) = Covariation(2,2) ;
    
    % + converging - pairwise
    [Corr,p,lb,ub] =  corrcoef(sn1_sim_res(:,a),circshift(sn2_sim_res(:,a),[minusPairwiseShift,0])) ; % correlations 
    sn_plusConverg_minusPairwise_lincoef(a) = Corr(1,2) ;
%    sn_plusConverg_minusPairwise_lincoefsig(a) = p(1,2) ;
%     if sn_plusConverg_minusPairwise_lincoefsig(a)>lincoefsigMin ;
%         sn_plusConverg_minusPairwise_lincoef(a) = 0 ;
%     end
%    sn_plusConverg_minusPairwise_std(a) = (Corr(1,2)-lb(1,2))*sqrt(length(sn1_sim_res(:,a)))/1.96 ; % std from confidence intervals
    Covariation =  cov(sn1_sim_res(:,a),circshift(sn2_sim_res(:,a),[minusPairwiseShift,0])) ; % covariation 
    sn_plusConverg_minusPairwise_cov(a) = Covariation(1,2) ;
    h=error_ellipse(Covariation,[sn1_sim_mean(a),sn2_sim_mean(a)],'conf',.95) ; % confidence elipse
    sn_plusConverg_minusPairwise_ellipseX(a,:) = get(h,'XData') ;
    sn_plusConverg_minusPairwise_ellipseY(a,:) = get(h,'YData') ;
    
    [Corr,p,lb,ub] =  corrcoef(sn1_shuff_res(:,a),circshift(sn2_shuff_res(:,a),[minusPairwiseShift,0])) ; % - converging - pairwise correlations 
    sn_minusConverg_minusPairwise_lincoef(a) = Corr(1,2) ;
%    sn_minusConverg_minusPairwise_lincoefsig(a) = p(1,2) ;
%     if sn_minusConverg_minusPairwise_lincoefsig(a)>lincoefsigMin ;
%         sn_minusConverg_minusPairwise_lincoef(a) = 0 ;
%     end
%    sn_minusConverg_minusPairwise_std(a) = (Corr(1,2)-lb(1,2))*sqrt(length(sn1_sim_res(:,a)))/1.96 ; % std from confidence intervals
    Covariation =  cov(sn1_shuff_res(:,a),circshift(sn2_shuff_res(:,a),[minusPairwiseShift,0])) ; % - converging - pairwise correlations 
    sn_minusConverg_minusPairwise_cov(a) = Covariation(1,2) ;    
    h=error_ellipse(Covariation,[sn1_shuff_mean(a),sn2_shuff_mean(a)],'conf',.95) ; % confidence elipse
    sn_minusConverg_minusPairwise_ellipseX(a,:) = get(h,'XData') ;
    sn_minusConverg_minusPairwise_ellipseY(a,:) = get(h,'YData') ;
    
    sn_corrCoefs(a,:) = [sn_plusConverg_plusPairwise_lincoef(a),sn_minusConverg_plusPairwise_lincoef(a),...
    sn_plusConverg_minusPairwise_lincoef(a),sn_minusConverg_minusPairwise_lincoef(a)] ;
    
%     sn_corrCoefsStd(a,:) = [sn_plusConverg_plusPairwise_std(a),sn_minusConverg_plusPairwise_std(a),...
%     sn_plusConverg_minusPairwise_std(a),sn_minusConverg_minusPairwise_std(a)] ;

    sn_covs(a,:) = [sn_plusConverg_plusPairwise_cov(a),sn_minusConverg_plusPairwise_cov(a),...
    sn_plusConverg_minusPairwise_cov(a),sn_minusConverg_minusPairwise_cov(a)] ;
end
close,clear h

% conductances
excg1_sim = nan(numTrials,epochSize) ;  % preallocate variables to improve memory usage and speed
excg2_sim = nan(numTrials,epochSize) ;   
excg1_shuff = nan(numTrials,epochSize) ;  
excg2_shuff = nan(numTrials,epochSize) ;  
inhg1_sim = nan(numTrials,epochSize) ; 
inhg2_sim = nan(numTrials,epochSize) ;  
inhg1_shuff = nan(numTrials,epochSize) ;  
inhg2_shuff = nan(numTrials,epochSize) ;  

for a = 1:numTrials ;
    excg1_sim(a,:) = excg(a*4-3,:) ; 
    excg2_sim(a,:) = excg(a*4-2,:) ;  
    excg1_shuff(a,:) = excg(a*4-1,:) ; 
    excg2_shuff(a,:) = excg(a*4,:) ; 
    
    inhg1_sim(a,:) = inhg(a*4-3,:) ; 
    inhg2_sim(a,:) = inhg(a*4-2,:) ; 
    inhg1_shuff(a,:) = inhg(a*4-1,:) ; 
    inhg2_shuff(a,:) = inhg(a*4,:) ; 
end

% mean g simultaneous
excg1_sim_mean = mean(excg1_sim) ;
inhg1_sim_mean = mean(inhg1_sim) ;
excg2_sim_mean = mean(excg2_sim) ;
inhg2_sim_mean = mean(inhg2_sim) ;

% variance g simultaneous
excg1_sim_var = var(excg1_sim) ;
inhg1_sim_var = var(inhg1_sim) ;
excg2_sim_var = var(excg2_sim) ;
inhg2_sim_var = var(inhg2_sim) ;

% shuffle check (are all the g waveforms the same in sim and shuff conditions?)
excg1_ShuffleCheck = min(excg1_sim_mean==mean(excg1_shuff)) ;
inhg1_ShuffleCheck = min(inhg1_sim_mean==mean(inhg1_shuff)) ;
excg2_ShuffleCheck = min(excg2_sim_mean==mean(excg2_shuff)) ;
inhg2_ShuffleCheck = min(inhg2_sim_mean==mean(inhg2_shuff)) ;
if min([excg1_ShuffleCheck,inhg1_ShuffleCheck,excg2_ShuffleCheck,inhg2_ShuffleCheck])<1;
   disp('not ideal shuffle')
end

% conductances divided by bar
for a=1:numBars ;
    excg1_sim_bar{a} = nan(numTrials,barPnts) ;  % preallocate variables to improve memory usage and speed 
    excg2_sim_bar{a} = nan(numTrials,barPnts) ;
    excg1_shuff_bar{a} = nan(numTrials,barPnts) ;
    excg2_shuff_bar{a} = nan(numTrials,barPnts) ;
    inhg1_sim_bar{a} = nan(numTrials,barPnts) ;
    inhg2_sim_bar{a} = nan(numTrials,barPnts) ;
    inhg1_shuff_bar{a} = nan(numTrials,barPnts) ;
    inhg2_shuff_bar{a} = nan(numTrials,barPnts) ;    
end

for a=1:numBars ;
    excg1_sim_bar{a} = excg1_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    excg2_sim_bar{a} = excg2_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    excg1_shuff_bar{a} = excg1_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    excg2_shuff_bar{a} = excg2_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    
    inhg1_sim_bar{a} = inhg1_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    inhg2_sim_bar{a} = inhg2_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    inhg1_shuff_bar{a} = inhg1_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    inhg2_shuff_bar{a} = inhg2_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;    
end

% figures

% figure % spike detection check
% for a=1:numTrials ; % on each trial
%     subplot(4,1,1)
%     plot(time,data(a*4-3,:))
%     hold on
%     plot(time,spikeTrain1_sim(a,:)*-30,'r')
%     hold off
%     
%     subplot(4,1,2)
%     plot(time,data(a*4-2,:))
%     hold on
%     plot(time,spikeTrain2_sim(a,:)*-30,'r')
%     hold off    
%         
%     subplot(4,1,3)
%     plot(time,data(a*4-1,:))
%     hold on
%     plot(time,spikeTrain1_shuff(a,:)*-30,'r')
%     hold off
%     
%     subplot(4,1,4)
%     plot(time,data(a*4,:))
%     hold on
%     plot(time,spikeTrain2_shuff(a,:)*-30,'r')
%     hold off    
%     
%     text(0,.9,num2str(epochs(a*4-3)),'units','norm')
%     pause
% end
% 
% pause
% close
% 
% figure % spike thresholds
% plot(SpikeThreshold_range,SpikeThreshold1_sim_Hist,'k')
% hold on
% plot(SpikeThreshold_range,SpikeThreshold2_sim_Hist,'k:')
% plot(SpikeThreshold_range,SpikeThreshold1_shuff_Hist,'g')
% hold on
% plot(SpikeThreshold_range,SpikeThreshold2_shuff_Hist,'g:')
% xlabel('spikeThreshold')
% ylabel('number of observations')
% legend('sim1','sim2','shuff1','shuff2')
% 
% pause
% close
% 
% 
% figure % injected conductances cell 1
% for a=1:numBars ;
%     subplot(2,numBars,a)
%     plot(time(1:barPnts),excg1_sim_bar{a},'b')
%     
%     subplot(2,numBars,a+numBars)
%     plot(time(1:barPnts),inhg1_sim_bar{a},'r')
% end
%     
% pause
% close
% 
% figure % injected conductances cell 2
% for a=1:numBars ;
%     subplot(2,numBars,a)
%     plot(time(1:barPnts),excg2_sim_bar{a},'b')
%     
%     subplot(2,numBars,a+numBars)
%     plot(time(1:barPnts),inhg2_sim_bar{a},'r')
% end
% 
% pause
% close
% 
% figure % injected conductances mean and variance
% subplot(2,1,1)
% plot(time,excg1_sim_mean,'b')
% hold on
% plot(time,excg1_sim_mean+excg1_sim_var,'b:')
% plot(time,excg1_sim_mean-excg1_sim_var,'b:')
% plot(time,inhg1_sim_mean,'r')
% plot(time,inhg1_sim_mean+inhg1_sim_var,'r:')
% plot(time,inhg1_sim_mean-inhg1_sim_var,'r:')
% 
% subplot(2,1,2)
% plot(time,excg2_sim_mean,'b')
% hold on
% plot(time,excg2_sim_mean+excg2_sim_var,'b:')
% plot(time,excg2_sim_mean-excg2_sim_var,'b:')
% plot(time,inhg2_sim_mean,'r')
% plot(time,inhg2_sim_mean+inhg2_sim_var,'r:')
% plot(time,inhg2_sim_mean-inhg2_sim_var,'r:')
% 
% 
% pause
% close
% 
figure % number of spikes per trial
for a = 1:numBars ;
    subplot(numBars,1,a)
    plot(sn1_sim(:,a),'k')
    hold on
    plot(sn2_sim(:,a),'k:')
    plot(sn1_shuff(:,a),'g')
    plot(sn2_shuff(:,a),'g:')
end
% 
% pause
% close
% 
% figure % spike number signal and noise correlations
% for a = 1:numBars ;
%     subplot(4,numBars,a)
%     plot(sn1_sim(:,a),sn2_sim(:,a),'k*')
%     
%     subplot(4,numBars,a+numBars)
%     plot(sn1_shuff(:,a),sn2_shuff(:,a),'g*')
%     
%     subplot(4,numBars,a+2*numBars)
%     plot(sn1_sim(:,a),circshift(sn2_sim(:,a),[2,0]),'ko')
%     
%     subplot(4,numBars,a+3*numBars)
%     plot(sn1_shuff(:,a),circshift(sn2_shuff(:,a),[2,0]),'go')
% end
% 
% pause
% close
% 
% figure % spike number noise correlations
% for a = 1:numBars ;
%     subplot(4,numBars,a)
%     plot(sn1_sim_res(:,a),sn2_sim_res(:,a),'k*')
%     
%     subplot(4,numBars,a+numBars)
%     plot(sn1_shuff_res(:,a),sn2_shuff_res(:,a),'g*')
%     
%     subplot(4,numBars,a+2*numBars)
%     plot(sn1_sim_res(:,a),circshift(sn2_sim_res(:,a),[2,0]),'ko')
%     
%     subplot(4,numBars,a+3*numBars)
%     plot(sn1_shuff_res(:,a),circshift(sn2_shuff_res(:,a),[2,0]),'go')
% end
% 
% pause
% close
% 
% 
% figure % mean and var number of spikes per bar
% subplot(2,1,1)
% errorbar([1:numBars],sn1_sim_mean,-sn1_sim_var,+sn1_sim_var,'k')
% hold on
% errorbar([1:numBars],sn1_shuff_mean,-sn1_shuff_var,+sn1_shuff_var,'g')
% title('cell 1')
% xlabel('bar')
% ylabel('spike number')
% 
% subplot(2,1,2)
% errorbar([1:numBars],sn2_sim_mean,-sn2_sim_var,+sn2_sim_var,'k')
% hold on
% errorbar([1:numBars],sn2_shuff_mean,-sn2_shuff_var,+sn2_shuff_var,'g')
% title('cell 2')
% xlabel('bar')
% ylabel('spike number')
% 
% pause
% close
% 
% figure % quantitative impact of converging and pairwise correlations on spike correlations
% for a = 1:numBars ;
%     subplot(1,numBars,a)
%     plot([1:4],sn_corrCoefs(a,:),'-o')
% 
%     set(gca,'XTick',[1:4],'XTickLabel',{'+c+p','-c+p','+c-p','-c-p'})
%     xlabel('input correlations')
%     ylabel('output sn correlation coef')
% end
% 
% pause
% close
% 
% figure % quantitative impact of converging and pairwise correlations on spike correlations
% for a = 1:numBars ;
%     subplot(1,numBars,a)
%     plot([1:4],sn_covs(a,:),'-o')
% 
%     set(gca,'XTick',[1:4],'XTickLabel',{'+c+p','-c+p','+c-p','-c-p'})
%     xlabel('input correlations')
%     ylabel('output sn covariance')
% end
% 
% pause
% close
% 
% 
% figure % bubble plots
% for a = 1:numBars ;
%     plot(sn_plusConverg_plusPairwise_ellipseX(a,:),sn_plusConverg_plusPairwise_ellipseY(a,:),'k')
%     hold on
%     plot(sn_minusConverg_plusPairwise_ellipseX(a,:),sn_minusConverg_plusPairwise_ellipseY(a,:),'g')
%     plot(sn_plusConverg_minusPairwise_ellipseX(a,:),sn_plusConverg_minusPairwise_ellipseY(a,:),'c')
%     plot(sn_minusConverg_minusPairwise_ellipseX(a,:),sn_minusConverg_minusPairwise_ellipseY(a,:),'y')
% end
% legend('+c+p','-c+p','+c-p','-c-p')

% FOR IGOR

% spike number and variance
% identifier = ['sn1SimMean',id,'cell',num2str(A)] ;
% ForIgor.(identifier) = sn1_sim_mean ; 
% 
% identifier = ['sn2SimMean',id,'cell',num2str(A)] ;
% ForIgor.(identifier) = sn2_sim_mean ; 
% 
% identifier = ['sn1SimVar',id,'cell',num2str(A)] ;
% ForIgor.(identifier) = sn1_sim_var ; 
% 
% identifier = ['sn2SimVar',id,'cell',num2str(A)] ;
% ForIgor.(identifier) = sn2_sim_var;

% spike number mean and variance per bar and condition
for a=1:numBars
    identifier = ['sn1SimMean',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn1_sim_mean(a) ; 

    identifier = ['sn2SimMean',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn2_sim_mean(a) ; 

    identifier = ['sn1SimVar',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn1_sim_var(a) ; 

    identifier = ['sn2SimVar',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn2_sim_var(a);

    identifier = ['sn1ShuffMean',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn1_shuff_mean(a) ; 

    identifier = ['sn2ShuffMean',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn2_shuff_mean(a) ; 

    identifier = ['sn1ShuffVar',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn1_shuff_var(a) ; 

    identifier = ['sn2ShuffVar',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn2_shuff_var(a);

end

% spike correlation coefficients and covariance by bar and condition
for a = 1:numBars ;
    identifier = ['snCorrCoefsPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_plusConverg_plusPairwise_lincoef(a);
    
    identifier = ['snCorrCoefsMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_minusConverg_plusPairwise_lincoef(a);
    
    identifier = ['snCorrCoefsMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_minusConverg_minusPairwise_lincoef(a);
    
    identifier = ['snCorrCoefsPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_plusConverg_minusPairwise_lincoef(a);
    
    identifier = ['snCovsPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_plusConverg_plusPairwise_cov(a);
    
    identifier = ['snCovsMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_minusConverg_plusPairwise_cov(a);
    
    identifier = ['snCovsPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_plusConverg_minusPairwise_cov(a);
    
    identifier = ['snCovsMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_minusConverg_minusPairwise_cov(a);
    
end



% make room in memmory
clearvars -except ForIgor Parameters id A excg* inhg* sT_corrCoefs sn_corrCoefs num* prePnts barPnts minusPairwiseShift SI time Voltage_mean SpikeThreshold_mean data SpikePnts absRef InhRev residualOption

%% spike generation and correlations through an LIF model
%{
samplerate = 1/SI ;
params.Eexc = 0 ;          % reversal potential for excitatory current
params.Einh = InhRev;
params.Eleak = Voltage_mean ;         % reversal potential of leak conductance (this can make cell spike spontaneously)

params.cap = 0.04 ; %nF                          % Set capacitance of cell
params.Vrest = Voltage_mean ;                 % initial potential of cell (resting potential mV) 
params.Vthresh = SpikeThreshold_mean ;     % spike threshold
params.AbsRef = absRef ; %sec       % absolute refractory period 

params.RelRefTau = .01 ; %sec              % decay time constant of relative refractory period
params.RelRefAmp = 4 ; %mV                  % amplitude of relative refractory period threshold change  

Iadd_range = [-500:10:100] ; % pA
Gleak_range = [0:10] ; % nS 
Vthresh_range = [SpikeThreshold_mean-10:SpikeThreshold_mean+10] ; % mV

% optimize parameters to match actual data
matlabpool 2 ; % open a pool of 2 matlab processes for parallel processesing
numberTrials_optimized = 60 ;
for b=1:numberTrials_optimized ;
    % optimize Iadd to capture the voltage mean accurately
    LIFv = nan(length(Iadd_range),length(excg(b,:))) ; % preallocate for speed
    parfor a=1:length(Iadd_range) ;
        LIFv(a,:) = LIFmodelGplusI(excg(b,:),inhg(b,:),0,Iadd_range(a),samplerate,params) ; % generate voltage data with LIF model
        MeanVerror(a) = (mean(LIFv(a,:)) - mean(data(b,:)))^2 ; 
    end
    [tempV,tempi] = min(MeanVerror) ;
    Iadd_opt(b) = Iadd_range(tempi) ;

    % optimize Gleak to capture the voltage variance accurately
    LIFv = nan(length(Gleak_range),length(excg(b,:))) ; % preallocate for speed
    parfor a=1:length(Gleak_range) ;
        LIFv(a,:) = LIFmodelGplusI(excg(b,:),inhg(b,:),Gleak_range(a),Iadd_opt(b),samplerate,params) ; % generate voltage data with LIF model
        varVerror(a) = (var(LIFv(a,:)) - var(data(b,:)))^2 ; 
    end
    [tempV,tempi] = min(varVerror) ;
    Gleak_opt(b) = Gleak_range(tempi) ;

    % optimize spike Thresh to capture the spike number accuratley
    LIFv = nan(length(Vthresh_range),length(excg(b,:))) ; % preallocate for speed
    for a =1:length(Vthresh_range) ;
        params.Vthresh = Vthresh_range(a) ;
        LIFv(a,:) = LIFmodelGplusI(excg(b,:),inhg(b,:),Gleak_opt(b),Iadd_opt(b),samplerate,params) ; 
        snerror(a) = (length(find(LIFv(a,:)==50)) - length(SpikePnts{b}))^2 ; 
    end
    [tempV,tempi] = min(snerror) ;
    Vthresh_opt(b) = Vthresh_range(tempi) ;
end
clear LIFv ;
matlabpool close

% spiking models
params.Vthresh = mean(Vthresh_opt) ;

% run LIF model with optimized parameters
LIFv_full = LIFmodelGplusI(excg,inhg,mean(Gleak_opt),mean(Iadd_opt),samplerate,params) ;

% run LIF model with spikes but linear synaptic integration
LIFv_reduced1 = LIFmodelGplusI_LinSynInt(excg,inhg,mean(Gleak_opt),mean(Iadd_opt),mean(Vthresh_opt),samplerate,params) ; 

% non spiking models
params.Vthresh = inf ;

% run LIF model with no spikes but standard nonlinear synaptic integration 
LIFv_noSpikes_reduced2 = LIFmodelGplusI(excg,inhg,mean(Gleak_opt),mean(Iadd_opt),samplerate,params) ;

% run LIF model with no spike threshold and linearize synaptic integration
LIFv_noSpikes_reduced3 = LIFmodelGplusI_LinSynInt(excg,inhg,mean(Gleak_opt),mean(Iadd_opt),mean(Vthresh_opt),samplerate,params) ;

% run LIF model with no spikes, linear synaptic integration, no Gleak or Iadd, and no voltage integration
LIFv_noSpikes_reduced4 = LIFmodelGplusI_LinSynInt(excg,inhg,0,0,mean(Vthresh_opt),samplerate,params) ;
LIFv_noSpikes_reduced4 = diff(LIFv_noSpikes_reduced4,[],2) ; % get rid of integration
LIFv_noSpikes_reduced4 = [LIFv_noSpikes_reduced4,LIFv_noSpikes_reduced4(:,end)] ;

for spikingModel = 1:2 ; % for each spiking model

    if spikingModel==1 ;
        LIFv = LIFv_full ;
    elseif spikingModel==2 ;
        LIFv = LIFv_reduced1 ;
    end

    % spike trains
    spikeTrain1_sim = zeros(numTrials,length(data)) ; 
    spikeTrain2_sim = zeros(numTrials,length(data)) ; 
    spikeTrain1_shuff = zeros(numTrials,length(data)) ;  
    spikeTrain2_shuff = zeros(numTrials,length(data)) ; 

    for a = 1:numTrials ;
        spikeTrain1_sim(a,:)=LIFv(a*4-3,:)==50 ; % mimic cell 1 + converging
        spikeTrain2_sim(a,:)=LIFv(a*4-2,:)==50 ; % mimic cell 2 + converging 
        spikeTrain1_shuff(a,:)=LIFv(a*4-1,:)==50 ; % mimic cell 1 - converging
        spikeTrain2_shuff(a,:)=LIFv(a*4,:)==50 ; % mimic cell 2 - converging
    end

    % spike trains divided by bar
    barPnts = floor((length(data))/numBars) ; 

    for a=1:numBars ; % preallocate variables to improve memory usage and speed
        spikeTrain1_sim_bar{a} = nan(numTrials,barPnts) ;
        spikeTrain2_sim_bar{a} = nan(numTrials,barPnts) ;
        spikeTrain1_shuff_bar{a} = nan(numTrials,barPnts) ;
        spikeTrain2_shuff_bar{a} = nan(numTrials,barPnts) ;
    end

    for a=1:numBars ;
        spikeTrain1_sim_bar{a} = spikeTrain1_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
        spikeTrain2_sim_bar{a} = spikeTrain2_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
        spikeTrain1_shuff_bar{a} = spikeTrain1_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
        spikeTrain2_shuff_bar{a} = spikeTrain2_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    end

    % spike number for each bar
    sn1_sim = nan(numTrials,numBars) ; % preallocate variables to improve memory usage and speed
    sn2_sim = nan(numTrials,numBars) ;
    sn1_shuff = nan(numTrials,numBars) ;
    sn2_shuff = nan(numTrials,numBars) ;
    sn1_sim_mean = nan(1,numBars) ;
    sn2_sim_mean = nan(1,numBars) ;
    sn1_shuff_mean = nan(1,numBars) ;
    sn2_shuff_mean = nan(1,numBars) ;

    for a=1:numBars ;
        sn1_sim(:,a) = sum(spikeTrain1_sim_bar{a},2) ;
        sn2_sim(:,a) = sum(spikeTrain2_sim_bar{a},2) ;
        sn1_shuff(:,a) = sum(spikeTrain1_shuff_bar{a},2) ;
        sn2_shuff(:,a) = sum(spikeTrain2_shuff_bar{a},2) ;

        sn1_sim_mean(a) = mean(sn1_sim(:,a)) ;
        sn2_sim_mean(a) = mean(sn2_sim(:,a)) ;
        sn1_shuff_mean(a) = mean(sn1_shuff(:,a)) ;
        sn2_shuff_mean(a) = mean(sn2_shuff(:,a)) ;
    end

    % spike number residuals
    sn1_sim_res = nan(numTrials-2*residualOption,numBars) ; % preallocate variables to improve memory usage and speed
    sn2_sim_res = nan(numTrials-2*residualOption,numBars) ;
    sn1_shuff_res = nan(numTrials-2*residualOption,numBars) ;
    sn2_shuff_res = nan(numTrials-2*residualOption,numBars) ;

    if residualOption == 0 ; % if full mean residual
        for a = 1:numTrials ;
            sn1_sim_res(a,:) = sn1_sim(a,:) - mean(sn1_sim) ;
            sn2_sim_res(a,:) = sn2_sim(a,:) - mean(sn2_sim) ;
            sn1_shuff_res(a,:) = sn1_shuff(a,:) - mean(sn1_shuff) ;
            sn2_shuff_res(a,:) = sn2_shuff(a,:) - mean(sn2_shuff) ;
        end
        
    elseif residualOption == 1 ; % if local residual
        for a = 2:numTrials-1 ;
            sn1_sim_res(a-1,:) = sn1_sim(a,:) - (sn1_sim(a-1,:)+sn1_sim(a+1,:))/2 ;
            sn2_sim_res(a-1,:) = sn2_sim(a,:) - (sn2_sim(a-1,:)+sn2_sim(a+1,:))/2 ;
            sn1_shuff_res(a-1,:) = sn1_shuff(a,:) - (sn1_shuff(a-1,:)+sn1_shuff(a+1,:))/2 ;
            sn2_shuff_res(a-1,:) = sn2_shuff(a,:) - (sn2_shuff(a-1,:)+sn2_shuff(a+1,:))/2 ;
        end
    end

    % spike number correlation (linear) for each bar
    lincoefsigMin = .05 ;
    for a=1:numBars ;
        [Corr,p] =  corrcoef(sn1_sim_res(:,a),sn2_sim_res(:,a)) ; % + converging + pairwise correlations 
        sn_plusConverg_plusPairwise_lincoef{spikingModel}(a) = Corr(1,2) ;
        sn_plusConverg_plusPairwise_lincoefsig(a) = p(1,2) ;
    %     if sn_plusConverg_plusPairwise_lincoefsig(a)>lincoefsigMin ;
    %         sn_plusConverg_plusPairwise_lincoef(a) = 0 ;
    %     end
        Covariation =  cov(sn1_sim_res(:,a),sn2_sim_res(:,a)) ; 
        sn_plusConverg_plusPairwise_cov{spikingModel}(a) = Covariation(1,2) ;
        sn1_sim_var{spikingModel}(a) = Covariation(1,1) ; % variance
        sn2_sim_var{spikingModel}(a) = Covariation(2,2) ;

        [Corr,p] =  corrcoef(sn1_shuff_res(:,a),sn2_shuff_res(:,a)) ; % - converging + pairwise correlations
        sn_minusConverg_plusPairwise_lincoef{spikingModel}(a) = Corr(1,2) ;
        sn_minusConverg_plusPairwise_lincoefsig(a) = p(1,2) ;
    %     if sn_minusConverg_plusPairwise_lincoefsig(a)>lincoefsigMin ;
    %         sn_minusConverg_plusPairwise_lincoef(a) = 0 ;
    %     end
        Covariation =  cov(sn1_shuff_res(:,a),sn2_shuff_res(:,a)) ; 
        sn_minusConverg_plusPairwise_cov{spikingModel}(a) = Covariation(1,2) ;
        sn1_shuff_var{spikingModel}(a) = Covariation(1,1) ; % variance
        sn2_shuff_var{spikingModel}(a) = Covariation(2,2) ;

        [Corr,p] =  corrcoef(sn1_sim_res(:,a),circshift(sn2_sim_res(:,a),[minusPairwiseShift,0])) ; % + converging - pairwise correlations 
        sn_plusConverg_minusPairwise_lincoef{spikingModel}(a) = Corr(1,2) ;
        sn_plusConverg_minusPairwise_lincoefsig(a) = p(1,2) ;
    %     if sn_plusConverg_minusPairwise_lincoefsig(a)>lincoefsigMin ;
    %         sn_plusConverg_minusPairwise_lincoef(a) = 0 ;
    %     end
        Covariation =  cov(sn1_sim_res(:,a),circshift(sn2_sim_res(:,a),[minusPairwiseShift,0])) ; 
        sn_plusConverg_minusPairwise_cov{spikingModel}(a) = Covariation(1,2) ;

        [Corr,p] =  corrcoef(sn1_shuff_res(:,a),circshift(sn2_shuff_res(:,a),[minusPairwiseShift,0])) ; % - converging - pairwise correlations 
        sn_minusConverg_minusPairwise_lincoef{spikingModel}(a) = Corr(1,2) ;
        sn_minusConverg_minusPairwise_lincoefsig(a) = p(1,2) ;
    %     if sn_minusConverg_minusPairwise_lincoefsig(a)>lincoefsigMin ;
    %         sn_minusConverg_minusPairwise_lincoef(a) = 0 ;
    %     end
        Covariation =  cov(sn1_shuff_res(:,a),circshift(sn2_shuff_res(:,a),[minusPairwiseShift,0])) ; 
        sn_minusConverg_minusPairwise_cov{spikingModel}(a) = Covariation(1,2) ;


        LIFsn_corrCoefs{spikingModel}(a,:) = [sn_plusConverg_plusPairwise_lincoef{spikingModel}(a),sn_minusConverg_plusPairwise_lincoef{spikingModel}(a),...
        sn_plusConverg_minusPairwise_lincoef{spikingModel}(a),sn_minusConverg_minusPairwise_lincoef{spikingModel}(a)] ;
    end

end

for nonSpikingModel = 1:3 ; % for each nonspiking model
    if nonSpikingModel == 1 ;
        LIFv_noSpikes = LIFv_noSpikes_reduced2 ;
    elseif nonSpikingModel == 2 ; 
        LIFv_noSpikes = LIFv_noSpikes_reduced3 ;
    elseif nonSpikingModel == 3 ;
        LIFv_noSpikes = LIFv_noSpikes_reduced4 ;
    end

    % voltages
    voltage1_sim = zeros(numTrials,length(LIFv_noSpikes)) ; 
    voltage2_sim = zeros(numTrials,length(LIFv_noSpikes)) ; 
    voltage1_shuff = zeros(numTrials,length(LIFv_noSpikes)) ;  
    voltage2_shuff = zeros(numTrials,length(LIFv_noSpikes)) ; 

    for a = 1:numTrials ;
        voltage1_sim(a,:)=LIFv_noSpikes(a*4-3,:) ; % mimic cell 1 + converging
        voltage2_sim(a,:)=LIFv_noSpikes(a*4-2,:) ; % mimic cell 2 + converging 
        voltage1_shuff(a,:)=LIFv_noSpikes(a*4-1,:) ; % mimic cell 1 - converging
        voltage2_shuff(a,:)=LIFv_noSpikes(a*4,:) ; % mimic cell 2 - converging
    end

    % voltage fluctuations divided by bar
    for a=1:numBars ; % preallocate variables to improve memory usage and speed
        voltage1_sim_bar{a} = nan(numTrials,barPnts) ;
        voltage2_sim_bar{a} = nan(numTrials,barPnts) ;
        voltage1_shuff_bar{a} = nan(numTrials,barPnts) ;
        voltage2_shuff_bar{a} = nan(numTrials,barPnts) ;
    end

    for a=1:numBars ;
        voltage1_sim_bar{a} = voltage1_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
        voltage2_sim_bar{a} = voltage2_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
        voltage1_shuff_bar{a} = voltage1_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
        voltage2_shuff_bar{a} = voltage2_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    end

    % mean voltage for each bar
    v1_sim = nan(numTrials,numBars) ; % preallocate variables to improve memory usage and speed
    v2_sim = nan(numTrials,numBars) ;
    v1_shuff = nan(numTrials,numBars) ;
    v2_shuff = nan(numTrials,numBars) ;

    for a=1:numBars ;
        v1_sim(:,a) = mean(voltage1_sim_bar{a},2) ;
        v2_sim(:,a) = mean(voltage2_sim_bar{a},2) ;
        v1_shuff(:,a) = mean(voltage1_shuff_bar{a},2) ;
        v2_shuff(:,a) = mean(voltage2_shuff_bar{a},2) ;
    end

    % mean voltage residuals
    v1_sim_res = nan(numTrials-2*residualOption,numBars) ; % preallocate variables to improve memory usage and speed
    v2_sim_res = nan(numTrials-2*residualOption,numBars) ;
    v1_shuff_res = nan(numTrials-2*residualOption,numBars) ;
    v2_shuff_res = nan(numTrials-2*residualOption,numBars) ;

    
    if residualOption == 0 ; % if full mean residual
            for a = 1:numTrials ;
                v1_sim_res(a,:) = v1_sim(a,:) - mean(v1_sim) ;
                v2_sim_res(a,:) = v2_sim(a,:) - mean(v2_sim) ;
                v1_shuff_res(a,:) = v1_shuff(a,:) - mean(v1_shuff) ;
                v2_shuff_res(a,:) = v2_shuff(a,:) - mean(v2_shuff) ;
            end
            
    elseif residualOption == 1 ; %local correction
        for a = 2:numTrials-1 ;
            v1_sim_res(a-1,:) = v1_sim(a,:) - (v1_sim(a-1,:)+v1_sim(a+1,:))/2 ;
            v2_sim_res(a-1,:) = v2_sim(a,:) - (v2_sim(a-1,:)+v2_sim(a+1,:))/2 ;
            v1_shuff_res(a-1,:) = v1_shuff(a,:) - (v1_shuff(a-1,:)+v1_shuff(a+1,:))/2 ;
            v2_shuff_res(a-1,:) = v2_shuff(a,:) - (v2_shuff(a-1,:)+v2_shuff(a+1,:))/2 ;
        end
    end


    % mean voltage correlation (linear) for each bar
    for a=1:numBars ;
        [Corr,p] =  corrcoef(v1_sim_res(:,a),v2_sim_res(:,a)) ; % + converging + pairwise correlations 
        v_plusConverg_plusPairwise_lincoef{nonSpikingModel}(a) = Corr(1,2) ;
        v_plusConverg_plusPairwise_lincoefsig(a) = p(1,2) ;
    %     if sn_plusConverg_plusPairwise_lincoefsig(a)>lincoefsigMin ;
    %         sn_plusConverg_plusPairwise_lincoef(a) = 0 ;
    %     end
        Covariation =  cov(v1_sim_res(:,a),v2_sim_res(:,a)) ; 
        v_plusConverg_plusPairwise_cov{nonSpikingModel}(a) = Covariation(1,2) ;
        v1_sim_var{nonSpikingModel}(a) = Covariation(1,1) ; % variance
        v2_sim_var{nonSpikingModel}(a) = Covariation(2,2) ;

        [Corr,p] =  corrcoef(v1_shuff_res(:,a),v2_shuff_res(:,a)) ; % - converging + pairwise correlations
        v_minusConverg_plusPairwise_lincoef{nonSpikingModel}(a) = Corr(1,2) ;
        v_minusConverg_plusPairwise_lincoefsig(a) = p(1,2) ;
    %     if sn_minusConverg_plusPairwise_lincoefsig(a)>lincoefsigMin ;
    %         sn_minusConverg_plusPairwise_lincoef(a) = 0 ;
    %     end
        Covariation =  cov(v1_shuff_res(:,a),v2_shuff_res(:,a)) ; 
        v_minusConverg_plusPairwise_cov{nonSpikingModel}(a) = Covariation(1,2) ;
        v1_shuff_var{nonSpikingModel}(a) = Covariation(1,1) ; % variance
        v2_shuff_var{nonSpikingModel}(a) = Covariation(2,2) ;

        [Corr,p] =  corrcoef(v1_sim_res(:,a),circshift(v2_sim_res(:,a),[minusPairwiseShift,0])) ; % + converging - pairwise correlations 
        v_plusConverg_minusPairwise_lincoef{nonSpikingModel}(a) = Corr(1,2) ;
        v_plusConverg_minusPairwise_lincoefsig(a) = p(1,2) ;
    %     if sn_plusConverg_minusPairwise_lincoefsig(a)>lincoefsigMin ;
    %         sn_plusConverg_minusPairwise_lincoef(a) = 0 ;
    %     end
        Covariation =  cov(v1_sim_res(:,a),circshift(v2_sim_res(:,a),[minusPairwiseShift,0])) ; 
        v_plusConverg_minusPairwise_cov{nonSpikingModel}(a) = Covariation(1,2) ;

        [Corr,p] =  corrcoef(v1_shuff_res(:,a),circshift(v2_shuff_res(:,a),[minusPairwiseShift,0])) ; % - converging - pairwise correlations 
        v_minusConverg_minusPairwise_lincoef{nonSpikingModel}(a) = Corr(1,2) ;
        v_minusConverg_minusPairwise_lincoefsig(a) = p(1,2) ;
    %     if sn_minusConverg_minusPairwise_lincoefsig(a)>lincoefsigMin ;
    %         sn_minusConverg_minusPairwise_lincoef(a) = 0 ;
    %     end
        Covariation =  cov(v1_shuff_res(:,a),circshift(v2_shuff_res(:,a),[minusPairwiseShift,0])) ; 
        v_minusConverg_minusPairwise_cov{nonSpikingModel}(a) = Covariation(1,2) ;

        LIFv_corrCoefs{nonSpikingModel}(a,:) = [v_plusConverg_plusPairwise_lincoef{nonSpikingModel}(a),v_minusConverg_plusPairwise_lincoef{nonSpikingModel}(a),...
        v_plusConverg_minusPairwise_lincoef{nonSpikingModel}(a),v_minusConverg_minusPairwise_lincoef{nonSpikingModel}(a)] ;
    end

end

% figure % comparing dynamic clamp and model voltage predictions
% for a=1:size(data,1) ;
%     plot(LIFv_full(a,:),'r')
%     hold on
%     plot(data(a,:))
%     hold off
%     pause
% end
% 
% pause
% close
% 
% figure % comparing dynamic clamp and model spike rasters 
% for a=1:size(data,1) ;
%     if ~isempty(find(SpikePnts{a},1))
%         plot(SpikePnts{a},a,'b*')
%     end
%     hold on
%     if ~isempty(find(LIFv_full(a,:)==50,1))
%         plot(find(LIFv_full(a,:)==50),a,'ro')
%     end
% end
% 
figure % comparing dynamic clamp and model spike number mean and variance
for a=1:length(SpikePnts) ;
    snData(a) = length(SpikePnts{a}) ;
    snLIF(a) = sum(LIFv_full(a,:)==50) ;
end
plot(snData,'b')
hold on
plot(snLIF,'r')


% FOR IGOR

% LIF for all conditions
% for a = 1:numBars ;
%     identifier = ['LIFsnCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = LIFsn_corrCoefs(a,:);
%     
%     identifier = ['LIFstCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = LIFsT_corrCoefs(a,:);    
% end

% % LIF for +c+p vs -c+p only
% numConditions = 2 ;
% 
% identifier = ['snCorrsActual',id,'cell',num2str(A)] ; % actual sn correlation
% ForIgor.(identifier) = sn_corrCoefs(1:numBars*numConditions) ;
%  
% identifier = ['snCorrsLIFfull',id,'cell',num2str(A)] ; % full lif model sn correlation
% ForIgor.(identifier) = LIFsn_corrCoefs{1}(1:numBars*numConditions);
%     
% identifier = ['snCorrsLIFreduced1',id,'cell',num2str(A)] ; % reduced lif model sn correlation
% ForIgor.(identifier) = LIFsn_corrCoefs{2}(1:numBars*numConditions);
% 
% identifier = ['vCorrsLIFreduced2',id,'cell',num2str(A)] ; % reduced lif model voltage correlation
% ForIgor.(identifier) = LIFv_corrCoefs{1}(1:numBars*numConditions);
% 
% identifier = ['vCorrsLIFreduced3',id,'cell',num2str(A)] ; % reduced lif model voltage correlation
% ForIgor.(identifier) = LIFv_corrCoefs{2}(1:numBars*numConditions);
% 
% identifier = ['vCorrsLIFreduced4',id,'cell',num2str(A)] ; % reduced lif model voltage correlation
% ForIgor.(identifier) = LIFv_corrCoefs{3}(1:numBars*numConditions);


% LIF spike correlation coefficients by bar and condition
for a = 1:numBars ;
    % full model
    identifier = ['LIFfsnCoefPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_plusConverg_plusPairwise_lincoef{1}(a);
    
    identifier = ['LIFfsnCoefMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_minusConverg_plusPairwise_lincoef{1}(a);
    
    identifier = ['LIFfsnCoefMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_minusConverg_minusPairwise_lincoef{1}(a);
    
    identifier = ['LIFfsnCoefPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_plusConverg_minusPairwise_lincoef{1}(a);
    
    % spiking reduced 1
    identifier = ['LIFr1snCoefPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_plusConverg_plusPairwise_lincoef{2}(a);
    
    identifier = ['LIFr1snCoefMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_minusConverg_plusPairwise_lincoef{2}(a);
    
    identifier = ['LIFr1snCoefMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_minusConverg_minusPairwise_lincoef{2}(a);
    
    identifier = ['LIFr1snCoefPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_plusConverg_minusPairwise_lincoef{2}(a);
    
    % nonspiking reduced 2
    identifier = ['LIFr2vCoefPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_plusConverg_plusPairwise_lincoef{1}(a);
    
    identifier = ['LIFr2vCoefMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_minusConverg_plusPairwise_lincoef{1}(a);
    
    identifier = ['LIFr2vCoefMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_minusConverg_minusPairwise_lincoef{1}(a);
    
    identifier = ['LIFr2vCoefPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_plusConverg_minusPairwise_lincoef{1}(a);
 
    % nonspiking reduced 3
    identifier = ['LIFr3vCoefPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_plusConverg_plusPairwise_lincoef{2}(a);
    
    identifier = ['LIFr3vCoefMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_minusConverg_plusPairwise_lincoef{2}(a);
    
    identifier = ['LIFr3vCoefMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_minusConverg_minusPairwise_lincoef{2}(a);
    
    identifier = ['LIFr3vCoefPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_plusConverg_minusPairwise_lincoef{2}(a);
    
    % nonspiking reduced 4
    identifier = ['LIFr4vCoefPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_plusConverg_plusPairwise_lincoef{3}(a);
    
    identifier = ['LIFr4vCoefMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_minusConverg_plusPairwise_lincoef{3}(a);
    
    identifier = ['LIFr4vCoefMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_minusConverg_minusPairwise_lincoef{3}(a);
    
    identifier = ['LIFr4vCoefPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_plusConverg_minusPairwise_lincoef{3}(a);

end

% LIF covariance by bar and condition
for a = 1:numBars ;
    % full model
    identifier = ['LIFfsnCovPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_plusConverg_plusPairwise_cov{1}(a);
    
    identifier = ['LIFfsnCovMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_minusConverg_plusPairwise_cov{1}(a);
    
    identifier = ['LIFfsnCovMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_minusConverg_minusPairwise_cov{1}(a);
    
    identifier = ['LIFfsnCovPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_plusConverg_minusPairwise_cov{1}(a);
    
    % spiking reduced 1
    identifier = ['LIFr1snCovPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_plusConverg_plusPairwise_cov{2}(a);
    
    identifier = ['LIFr1snCovMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_minusConverg_plusPairwise_cov{2}(a);
    
    identifier = ['LIFr1snCovMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_minusConverg_minusPairwise_cov{2}(a);
    
    identifier = ['LIFr1snCovPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_plusConverg_minusPairwise_cov{2}(a);
    
    % nonspiking reduced 2
    identifier = ['LIFr2vCovPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_plusConverg_plusPairwise_cov{1}(a);
    
    identifier = ['LIFr2vCovMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_minusConverg_plusPairwise_cov{1}(a);
    
    identifier = ['LIFr2vCovMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_minusConverg_minusPairwise_cov{1}(a);
    
    identifier = ['LIFr2vCovPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_plusConverg_minusPairwise_cov{1}(a);
 
    % nonspiking reduced 3
    identifier = ['LIFr3vCovPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_plusConverg_plusPairwise_cov{2}(a);
    
    identifier = ['LIFr3vCovMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_minusConverg_plusPairwise_cov{2}(a);
    
    identifier = ['LIFr3vCovMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_minusConverg_minusPairwise_cov{2}(a);
    
    identifier = ['LIFr3vCovPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_plusConverg_minusPairwise_cov{2}(a);
    
    % nonspiking reduced 4
    identifier = ['LIFr4vCovPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_plusConverg_plusPairwise_cov{3}(a);
    
    identifier = ['LIFr4vCovMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_minusConverg_plusPairwise_cov{3}(a);
    
    identifier = ['LIFr4vCovMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_minusConverg_minusPairwise_cov{3}(a);
    
    identifier = ['LIFr4vCovPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = v_plusConverg_minusPairwise_cov{3}(a);

end

% LIF variance by bar and condition
for a = 1:numBars ;
    % full model
    identifier = ['LIFfsnVar1Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  sn1_sim_var{1}(a) ;
    
    identifier = ['LIFfsnVar2Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  sn2_sim_var{1}(a) ;
 
    identifier = ['LIFfsnVar1Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  sn1_shuff_var{1}(a) ;
    
    identifier = ['LIFfsnVar2Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  sn2_shuff_var{1}(a) ;
    
    % spiking reduced 1
    identifier = ['LIFr1snVar1Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  sn1_sim_var{2}(a) ;
    
    identifier = ['LIFr1snVar2Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  sn2_sim_var{2}(a) ;
 
    identifier = ['LIFr1snVar1Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  sn1_shuff_var{2}(a) ;
    
    identifier = ['LIFr1snVar2Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  sn2_shuff_var{2}(a) ;
    
    % nonspiking reduced 2
    identifier = ['LIFr2vVar1Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  v1_sim_var{1}(a) ;
    
    identifier = ['LIFr2vVar2Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  v2_sim_var{1}(a) ;
 
    identifier = ['LIFr2vVar1Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  v1_shuff_var{1}(a) ;
    
    identifier = ['LIFr2vVar2Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  v2_shuff_var{1}(a) ;
 
    % nonspiking reduced 3
    identifier = ['LIFr3vVar1Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  v1_sim_var{2}(a) ;
    
    identifier = ['LIFr3vVar2Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  v2_sim_var{2}(a) ;
 
    identifier = ['LIFr3vVar1Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  v1_shuff_var{2}(a) ;
    
    identifier = ['LIFr3vVar2Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  v2_shuff_var{2}(a) ;
    
    % nonspiking reduced 4
    identifier = ['LIFr4vVar1Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  v1_sim_var{3}(a) ;
    
    identifier = ['LIFr4vVar2Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  v2_sim_var{3}(a) ;
 
    identifier = ['LIFr4vVar1Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  v1_shuff_var{3}(a) ;
    
    identifier = ['LIFr4vVar2Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) =  v2_shuff_var{3}(a) ;
end

%}
% make room in memmory
clearvars -except ForIgor Parameters id A excg* inhg* sT_corrCoefs sn_corrCoefs num* prePnts barPnts minusPairwiseShift SI time Voltage_mean InhRev SpikeThreshold_mean residualOption Vthresh_opt

%% noise correlations and variance predicted from the presented conductance means  
%{
% each conductance trial averaged over time
for b=1:numBars ;
    excg1_sim_bar_mean{b} = mean(excg1_sim_bar{b},2) ;
    excg2_sim_bar_mean{b} = mean(excg2_sim_bar{b},2) ;
    excg1_shuff_bar_mean{b} = mean(excg1_shuff_bar{b},2) ;
    excg2_shuff_bar_mean{b} = mean(excg2_shuff_bar{b},2) ;

    inhg1_sim_bar_mean{b} = mean(inhg1_sim_bar{b},2) ;
    inhg2_sim_bar_mean{b} = mean(inhg2_sim_bar{b},2) ;
    inhg1_shuff_bar_mean{b} = mean(inhg1_shuff_bar{b},2) ;
    inhg2_shuff_bar_mean{b} = mean(inhg2_shuff_bar{b},2) ;
end

% conductance residuals
for b=1:numBars ; % preallocate variables to improve memory usage and speed 
    excg1_sim_bar_res{b} = nan(numTrials-2*residualOption,1) ;
    excg2_sim_bar_res{b} = nan(numTrials-2*residualOption,1) ;
    excg1_shuff_bar_res{b} = nan(numTrials-2*residualOption,1) ;
    excg2_shuff_bar_res{b} = nan(numTrials-2*residualOption,1) ;
    inhg1_sim_bar_res{b} = nan(numTrials-2*residualOption,1) ;
    inhg2_sim_bar_res{b} = nan(numTrials-2*residualOption,1) ;
    inhg1_shuff_bar_res{b} = nan(numTrials-2*residualOption,1) ;
    inhg2_shuff_bar_res{b} = nan(numTrials-2*residualOption,1) ;
end

for b=1:numBars ;
    if residualOption == 0 ; % if full mean residual
        for a = 1:numTrials ; % non local
            excg1_sim_bar_res{b}(a,:) = excg1_sim_bar_mean{b}(a) - mean(excg1_sim_bar_mean{b}) ;
            excg2_sim_bar_res{b}(a,:) = excg2_sim_bar_mean{b}(a) - mean(excg2_sim_bar_mean{b}) ;
            excg1_shuff_bar_res{b}(a,:) = excg1_shuff_bar_mean{b}(a) - mean(excg1_shuff_bar_mean{b}) ;
            excg2_shuff_bar_res{b}(a,:) = excg2_shuff_bar_mean{b}(a) - mean(excg2_shuff_bar_mean{b}) ;

            inhg1_sim_bar_res{b}(a,:) = inhg1_sim_bar_mean{b}(a) - mean(inhg1_sim_bar_mean{b}) ;
            inhg2_sim_bar_res{b}(a,:) = inhg2_sim_bar_mean{b}(a) - mean(inhg2_sim_bar_mean{b}) ;
            inhg1_shuff_bar_res{b}(a,:) = inhg1_shuff_bar_mean{b}(a) - mean(inhg1_shuff_bar_mean{b}) ;
            inhg2_shuff_bar_res{b}(a,:) = inhg2_shuff_bar_mean{b}(a) - mean(inhg2_shuff_bar_mean{b}) ;
        end
    elseif residualOption == 1 ; % if local correction
        for a = 2:numTrials-1 ;  % 
            excg1_sim_bar_res{b}(a-1,:) = excg1_sim_bar_mean{b}(a) - (excg1_sim_bar_mean{b}(a-1)+excg1_sim_bar_mean{b}(a+1))/2 ;
            excg2_sim_bar_res{b}(a-1,:) = excg2_sim_bar_mean{b}(a) - (excg2_sim_bar_mean{b}(a-1)+excg2_sim_bar_mean{b}(a+1))/2 ;
            excg1_shuff_bar_res{b}(a-1,:) = excg1_shuff_bar_mean{b}(a) - (excg1_shuff_bar_mean{b}(a-1)+excg1_shuff_bar_mean{b}(a+1))/2 ;
            excg2_shuff_bar_res{b}(a-1,:) = excg2_shuff_bar_mean{b}(a) - (excg2_shuff_bar_mean{b}(a-1)+excg2_shuff_bar_mean{b}(a+1))/2 ;

            inhg1_sim_bar_res{b}(a-1,:) = inhg1_sim_bar_mean{b}(a) - (inhg1_sim_bar_mean{b}(a-1)+inhg1_sim_bar_mean{b}(a+1))/2 ;
            inhg2_sim_bar_res{b}(a-1,:) = inhg2_sim_bar_mean{b}(a) - (inhg2_sim_bar_mean{b}(a-1)+inhg2_sim_bar_mean{b}(a+1))/2 ;
            inhg1_shuff_bar_res{b}(a-1,:) = inhg1_shuff_bar_mean{b}(a) - (inhg1_shuff_bar_mean{b}(a-1)+inhg1_shuff_bar_mean{b}(a+1))/2 ;
            inhg2_shuff_bar_res{b}(a-1,:) = inhg2_shuff_bar_mean{b}(a) - (inhg2_shuff_bar_mean{b}(a-1)+inhg2_shuff_bar_mean{b}(a+1))/2 ;
        end  
    end
end

% g correlations
for b=1:numBars ;
    % +ei +p
    temp = cov(excg1_sim_bar_res{b}, inhg1_sim_bar_res{b}) ;
    excg1_sim_bar_var(b) = temp(1,1) ;
    inhg1_sim_bar_var(b) = temp(2,2) ;
    excg1inhg1_sim_bar_cov(b) = temp(1,2) ;
    
    temp = cov(excg2_sim_bar_res{b}, inhg2_sim_bar_res{b}) ;
    excg2_sim_bar_var(b) = temp(1,1) ;
    inhg2_sim_bar_var(b) = temp(2,2) ;
    excg2inhg2_sim_bar_cov(b) = temp(1,2) ;   
    
    temp = cov(excg1_sim_bar_res{b}, excg2_sim_bar_res{b}) ;
    excg1excg2_sim_bar_cov(b) = temp(1,2) ;   
    
    temp = cov(inhg1_sim_bar_res{b}, inhg2_sim_bar_res{b}) ;
    inhg1inhg2_sim_bar_cov(b) = temp(1,2) ;    
    
    temp = cov(excg1_sim_bar_res{b}, inhg2_sim_bar_res{b}) ;
    excg1inhg2_sim_bar_cov(b) = temp(1,2) ;   
    
    temp = cov(inhg1_sim_bar_res{b}, excg2_sim_bar_res{b}) ;
    inhg1excg2_sim_bar_cov(b) = temp(1,2) ; 
    
    % -ei +p
    temp = cov(excg1_shuff_bar_res{b}, inhg1_shuff_bar_res{b}) ;
    excg1_shuff_bar_var(b) = temp(1,1) ;
    inhg1_shuff_bar_var(b) = temp(2,2) ;
    excg1inhg1_shuff_bar_cov(b) = temp(1,2) ;
    
    temp = cov(excg2_shuff_bar_res{b}, inhg2_shuff_bar_res{b}) ;
    excg2_shuff_bar_var(b) = temp(1,1) ;
    inhg2_shuff_bar_var(b) = temp(2,2) ;
    excg2inhg2_shuff_bar_cov(b) = temp(1,2) ;
    
    temp = cov(excg1_shuff_bar_res{b}, excg2_shuff_bar_res{b}) ;
    excg1excg2_shuff_bar_cov(b) = temp(1,2) ;   
    
    temp = cov(inhg1_shuff_bar_res{b}, inhg2_shuff_bar_res{b}) ;
    inhg1inhg2_shuff_bar_cov(b) = temp(1,2) ; 
    
    temp = cov(excg1_shuff_bar_res{b}, inhg2_shuff_bar_res{b}) ;
    excg1inhg2_shuff_bar_cov(b) = temp(1,2) ;   
    
    temp = cov(inhg1_shuff_bar_res{b}, excg2_shuff_bar_res{b}) ;
    inhg1excg2_shuff_bar_cov(b) = temp(1,2) ; 
    
    % +ei -p
    temp = cov(excg1_sim_bar_res{b}, circshift(excg2_sim_bar_res{b},minusPairwiseShift)) ;
    excg1excg2Shift_sim_bar_cov(b) = temp(1,2) ;   
    
    temp = cov(inhg1_sim_bar_res{b}, circshift(inhg2_sim_bar_res{b},minusPairwiseShift)) ;
    inhg1inhg2Shift_sim_bar_cov(b) = temp(1,2) ;  
    
    temp = cov(excg1_sim_bar_res{b}, circshift(inhg2_sim_bar_res{b},minusPairwiseShift)) ;
    excg1inhg2Shift_sim_bar_cov(b) = temp(1,2) ;   
    
    temp = cov(inhg1_sim_bar_res{b},circshift(excg2_sim_bar_res{b},minusPairwiseShift)) ;
    inhg1excg2Shift_sim_bar_cov(b) = temp(1,2) ; 
    
    % -ei -p
    temp = cov(excg1_shuff_bar_res{b}, circshift(excg2_shuff_bar_res{b},minusPairwiseShift)) ;
    excg1excg2Shift_shuff_bar_cov(b) = temp(1,2) ;   
    
    temp = cov(inhg1_shuff_bar_res{b}, circshift(inhg2_shuff_bar_res{b},minusPairwiseShift)) ;
    inhg1inhg2Shift_shuff_bar_cov(b) = temp(1,2) ;  
    
    temp = cov(excg1_shuff_bar_res{b}, circshift(inhg2_shuff_bar_res{b},minusPairwiseShift)) ;
    excg1inhg2Shift_shuff_bar_cov(b) = temp(1,2) ;   
    
    temp = cov(inhg1_shuff_bar_res{b},circshift(excg2_shuff_bar_res{b},minusPairwiseShift)) ;
    inhg1excg2Shift_shuff_bar_cov(b) = temp(1,2) ; 
end

% linear analytic prediction

% analytical solution assuming exc - inh is proportional to spike output  
% Pss = (Ve1e2 +a^2*Vi1i2 -a*Ve1i2 -a*Ve2i1)/sqrt((Ve1^2 +a^2Vi1^2 -2aVe1i1)*(Ve2^2 +a^2Vi2^2 -2aVe2i2))
% a (alpha) is ratio of inh/exc
Ve1e2 = [excg1excg2_sim_bar_cov',excg1excg2_shuff_bar_cov',excg1excg2Shift_sim_bar_cov',excg1excg2Shift_shuff_bar_cov'] ;
Vi1i2 = [inhg1inhg2_sim_bar_cov',inhg1inhg2_shuff_bar_cov',inhg1inhg2Shift_sim_bar_cov',inhg1inhg2Shift_shuff_bar_cov'] ;
Ve1i2 = [excg1inhg2_sim_bar_cov',excg1inhg2_shuff_bar_cov',excg1inhg2Shift_sim_bar_cov',excg1inhg2Shift_shuff_bar_cov'] ;
Vi1e2 = [inhg1excg2_sim_bar_cov',inhg1excg2_shuff_bar_cov',inhg1excg2Shift_sim_bar_cov',inhg1excg2Shift_shuff_bar_cov'] ;

Ve1 = [excg1_sim_bar_var',excg1_shuff_bar_var',excg1_sim_bar_var',excg1_shuff_bar_var'] ;
Ve2 = [excg2_sim_bar_var',excg2_shuff_bar_var',excg2_sim_bar_var',excg2_shuff_bar_var'] ;
Vi1 = [inhg1_sim_bar_var',inhg1_shuff_bar_var',inhg1_sim_bar_var',inhg1_shuff_bar_var'] ;
Vi2 = [inhg2_sim_bar_var',inhg2_shuff_bar_var',inhg2_sim_bar_var',inhg2_shuff_bar_var'] ;
Ve1i1 = [excg1inhg1_sim_bar_cov',excg1inhg1_shuff_bar_cov',excg1inhg1_sim_bar_cov',excg1inhg1_shuff_bar_cov'] ;
Ve2i2 = [excg2inhg2_sim_bar_cov',excg2inhg2_shuff_bar_cov',excg2inhg2_sim_bar_cov',excg2inhg2_shuff_bar_cov'] ;

numConditions = 4 ; 
alpha =  abs((mean(Vthresh_opt)-InhRev)/(mean(Vthresh_opt)-0)) ; % alpha from LIF spike threshold
for a=1:numBars; % for every bar
    for b=1:numConditions ; % for every set pairings   
        sn_covariance_predicted(a,b) = Ve1e2(a,b) + alpha^2*Vi1i2(a,b) - alpha*Ve1i2(a,b) - alpha*Vi1e2(a,b) ;
        sn_variance1_predicted(a,b) = Ve1(a,b) + alpha^2*Vi1(a,b) - 2*alpha*Ve1i1(a,b) ;
        sn_variance2_predicted(a,b) = Ve2(a,b) + alpha^2*Vi2(a,b) - 2*alpha*Ve2i2(a,b) ;
        sn_corrcoef_predicted(a,b) = sn_covariance_predicted(a,b)/sqrt(sn_variance1_predicted(a,b)*sn_variance2_predicted(a,b)) ;
    end            
end





% figures
%{
figure % estimated alpha values
plot(alphaHistX,alpha_snHist,'o-')
hold on
plot(alphaHistX,alpha_stHist,'+-')
xlabel('alpha')
ylabel('number of observations')


figure % actual vs predicted
subplot(1,3,1)
plot(alpha_range,CorrActPred_sn,'o-')
hold on
plot(alpha_range,CorrActPredNoCross_sn,'o--')
plot(alpha_range,CorrActPred_st,'+-')
hold on
plot(alpha_range,CorrActPredNoCross_st,'+--')
xlabel('alpha')
ylabel('corr act v pred')

subplot(1,3,2)
plot(alpha_range,Mse_sn,'o-')
hold on
plot(alpha_range,MseNoCross_sn,'o--')
plot(alpha_range,Mse_st,'+-')
hold on
plot(alpha_range,MseNoCross_st,'+--')
xlabel('alpha')
ylabel('mse')

subplot(1,3,3)
plot(sn_corrCoefs(1:numel(Pss_predicted_sn)),Pss_predicted_sn(:),'*')
hold on
xlabel('actual')
ylabel('predicted')
legend('sn','st')

%}

% FOR IGOR
%
% identifier = ['Ve1e2',id,'cell',num2str(A)] ; % Ve1e2 
% ForIgor.(identifier) = Ve1e2 ; 
% 
% identifier = ['Vi1i2',id,'cell',num2str(A)] ; % Vi1i2 
% ForIgor.(identifier) = Vi1i2 ; 
% 
% identifier = ['Ve1i2',id,'cell',num2str(A)] ; % Ve1i2 
% ForIgor.(identifier) = Ve1i2 ; 
% 
% identifier = ['Vi1e2',id,'cell',num2str(A)] ; % Vi1e2 
% ForIgor.(identifier) = Vi1e2 ;
% 
% identifier = ['Ve1i1',id,'cell',num2str(A)] ; % Ve1i1
% ForIgor.(identifier) = Ve1i1 ; 
% 
% identifier = ['Ve2i2',id,'cell',num2str(A)] ; % Ve2i2
% ForIgor.(identifier) = Ve2i2 ; 
% 
% identifier = ['Ve1',id,'cell',num2str(A)] ; % Ve1
% ForIgor.(identifier) = Ve1 ; 
% 
% identifier = ['Vi1',id,'cell',num2str(A)] ; % Vi1
% ForIgor.(identifier) = Vi1 ; 
% 
% identifier = ['Ve2',id,'cell',num2str(A)] ; % Ve2
% ForIgor.(identifier) = Ve1 ; 
% 
% identifier = ['Vi2',id,'cell',num2str(A)] ; % Vi2
% ForIgor.(identifier) = Vi2 ; 

%

% identifier = ['snCorrsActual',id,'cell',num2str(A)] ;
% ForIgor.(identifier) = sn_corrCoefs(1:numel(Pss_predicted_sn)) ;

% identifier = ['snCorrsAnalytic',id,'cell',num2str(A)] ;
% ForIgor.(identifier) = Pss_predicted_sn(:) ;
%
% identifier = ['sTCorrsActual',id,'cell',num2str(A)] ;
% ForIgor.(identifier) = sT_corrCoefs(1:numel(Pss_predicted)) ;
% 
% identifier = ['sTCorrsAnalytic',id,'cell',num2str(A)] ;
% ForIgor.(identifier) = Pss_predicted_st(:) ;
% 
% for a = 1:numBars ;
%     identifier = ['linAnalyticsnCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = Pss_predicted_sn(a,:) ;     
% end
    
% identifier = ['eMinusIcorr',id,'cell',num2str(A)] ;
% ForIgor.(identifier) = eMinusI_corrCoefs(1:numel(Pss_predicted_sn)) ;
    
% linear anlaytic predicted spike correlation coefficients by bar and condition
for a = 1:numBars ;
    identifier = ['LapCoefsPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_corrcoef_predicted(a,1);
    
    identifier = ['LapCoefsMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_corrcoef_predicted(a,2);

    identifier = ['LapCoefsPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_corrcoef_predicted(a,3);
    
    identifier = ['LapCoefsMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_corrcoef_predicted(a,4);
    
    
    identifier = ['LapCovsPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_covariance_predicted(a,1);
    
    identifier = ['LapCovsMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_covariance_predicted(a,2);

    identifier = ['LapCovsPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_covariance_predicted(a,3);
    
    identifier = ['LapCovsMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_covariance_predicted(a,4);
    
 
    
    identifier = ['LapVar1Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_variance1_predicted(a,1);
    
    identifier = ['LapVar1Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_variance1_predicted(a,2);
    
    identifier = ['LapVar2Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_variance2_predicted(a,1);
    
    identifier = ['LapVar2Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_variance2_predicted(a,2);
end

% % input covariance for each bar and condition
% for a = 1:numBars ;
%     % +ei +p
%     identifier = ['CovE1E2PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1excg2_sim_bar_cov(a) ;
% 
%     identifier = ['CovI1I2PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1inhg2_sim_bar_cov(a) ;
%     
%     identifier = ['CovE1I2PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1inhg2_sim_bar_cov(a) ;    
% 
%     identifier = ['CovE2I1PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1excg2_sim_bar_cov(a) ;
%    
%     identifier = ['CovE1I1PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1inhg1_sim_bar_cov(a) ;
%     
%     identifier = ['CovE2I2PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg2inhg2_sim_bar_cov(a) ;
%     
%     % -ei +p
%     identifier = ['CovE1E2McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1excg2_shuff_bar_cov(a) ;
% 
%     identifier = ['CovI1I2McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1inhg2_shuff_bar_cov(a) ;
%     
%     identifier = ['CovE1I2McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1inhg2_shuff_bar_cov(a) ;    
% 
%     identifier = ['CovE2I1McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1excg2_shuff_bar_cov(a) ;
%    
%     identifier = ['CovE1I1McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1inhg1_shuff_bar_cov(a) ;
%     
%     identifier = ['CovE2I2McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg2inhg2_shuff_bar_cov(a) ;
%     
%     % -ei -p
%     identifier = ['CovE1E2McMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1excg2Shift_shuff_bar_cov(a) ;
% 
%     identifier = ['CovI1I2McMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1inhg2Shift_shuff_bar_cov(a) ;
%     
%     identifier = ['CovE1I2McMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1inhg2Shift_shuff_bar_cov(a) ;    
% 
%     identifier = ['CovE2I1McMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1excg2Shift_shuff_bar_cov(a) ;
%    
%     % +ei -p
%     identifier = ['CovE1E2PcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1excg2Shift_sim_bar_cov(a) ;
% 
%     identifier = ['CovI1I2PcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1inhg2Shift_sim_bar_cov(a) ;
%     
%     identifier = ['CovE1I2PcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1inhg2Shift_sim_bar_cov(a) ;    
% 
%     identifier = ['CovE2I1PcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1excg2Shift_sim_bar_cov(a) ;
% 
% end

%}
% make room in memmory
clearvars -except ForIgor Parameters id A excg* inhg* sn_corrCoefs num* prePnts barPnts minusPairwiseShift SI time* residualOption Vthresh_opt InhRev


%% noise correlations and variance predicted from the presented conductance cross corr peaks  
%{
Vthresh_opt=-40 ; % TEMP WHEN NOT RUNNING LIF!!!!!!!!
% conductance residuals (local correction)
for b=1:numBars ; % preallocate variables to improve memory usage and speed 
    excg1_sim_bar_res{b} = nan(numTrials-2*residualOption,barPnts) ;
    excg2_sim_bar_res{b} = nan(numTrials-2*residualOption,barPnts) ;
    excg1_shuff_bar_res{b} = nan(numTrials-2*residualOption,barPnts) ;
    excg2_shuff_bar_res{b} = nan(numTrials-2*residualOption,barPnts) ;
    inhg1_sim_bar_res{b} = nan(numTrials-2*residualOption,barPnts) ;
    inhg2_sim_bar_res{b} = nan(numTrials-2*residualOption,barPnts) ;
    inhg1_shuff_bar_res{b} = nan(numTrials-2*residualOption,barPnts) ;
    inhg2_shuff_bar_res{b} = nan(numTrials-2*residualOption,barPnts) ;
end

for b=1:numBars ;
    if residualOption==0 ;
        for a = 1:numTrials ;
            excg1_sim_bar_res{b}(a,:) = excg1_sim_bar{b}(a,:) - mean(excg1_sim_bar_res{b}) ;
            excg2_sim_bar_res{b}(a,:) = excg2_sim_bar{b}(a,:) - mean(excg2_sim_bar_res{b}) ;
            excg1_shuff_bar_res{b}(a,:) = excg1_shuff_bar{b}(a,:) - mean(excg1_shuff_bar_res{b}) ;
            excg2_shuff_bar_res{b}(a,:) = excg2_shuff_bar{b}(a,:) - mean(excg2_shuff_bar_res{b}) ;

            inhg1_sim_bar_res{b}(a,:) = inhg1_sim_bar{b}(a,:) - mean(inhg1_sim_bar_res{b}) ;
            inhg2_sim_bar_res{b}(a,:) = inhg2_sim_bar{b}(a,:) - mean(inhg2_sim_bar_res{b}) ;
            inhg1_shuff_bar_res{b}(a,:) = inhg1_shuff_bar{b}(a,:) - mean(inhg1_shuff_bar_res{b}) ;
            inhg2_shuff_bar_res{b}(a,:) = inhg2_shuff_bar{b}(a,:) - mean(inhg2_shuff_bar_res{b}) ;
        end
        
    elseif residualOption==1 ;
        for a = 2:numTrials-1 ;
            excg1_sim_bar_res{b}(a-1,:) = excg1_sim_bar{b}(a,:) - (excg1_sim_bar{b}(a-1,:)+excg1_sim_bar{b}(a+1,:))/2 ;
            excg2_sim_bar_res{b}(a-1,:) = excg2_sim_bar{b}(a,:) - (excg2_sim_bar{b}(a-1,:)+excg2_sim_bar{b}(a+1,:))/2 ;
            excg1_shuff_bar_res{b}(a-1,:) = excg1_shuff_bar{b}(a,:) - (excg1_shuff_bar{b}(a-1,:)+excg1_shuff_bar{b}(a+1,:))/2 ;
            excg2_shuff_bar_res{b}(a-1,:) = excg2_shuff_bar{b}(a,:) - (excg2_shuff_bar{b}(a-1,:)+excg2_shuff_bar{b}(a+1,:))/2 ;

            inhg1_sim_bar_res{b}(a-1,:) = inhg1_sim_bar{b}(a,:) - (inhg1_sim_bar{b}(a-1,:)+inhg1_sim_bar{b}(a+1,:))/2 ;
            inhg2_sim_bar_res{b}(a-1,:) = inhg2_sim_bar{b}(a,:) - (inhg2_sim_bar{b}(a-1,:)+inhg2_sim_bar{b}(a+1,:))/2 ;
            inhg1_shuff_bar_res{b}(a-1,:) = inhg1_shuff_bar{b}(a,:) - (inhg1_shuff_bar{b}(a-1,:)+inhg1_shuff_bar{b}(a+1,:))/2 ;
            inhg2_shuff_bar_res{b}(a-1,:) = inhg2_shuff_bar{b}(a,:) - (inhg2_shuff_bar{b}(a-1,:)+inhg2_shuff_bar{b}(a+1,:))/2 ;
        end
    end
    
    % biased variance 
    excg1_sim_bar_var(b) = mean(mean(excg1_sim_bar_res{b}.^2,2)) ; 
    excg2_sim_bar_var(b) = mean(mean(excg2_sim_bar_res{b}.^2,2)) ;
    excg1_shuff_bar_var(b) = mean(mean(excg1_shuff_bar_res{b}.^2,2)) ;
    excg2_shuff_bar_var(b) = mean(mean(excg2_shuff_bar_res{b}.^2,2)) ;

    inhg1_sim_bar_var(b) = mean(mean(inhg1_sim_bar_res{b}.^2,2)) ;
    inhg2_sim_bar_var(b) = mean(mean(inhg2_sim_bar_res{b}.^2,2)) ;
    inhg1_shuff_bar_var(b) = mean(mean(inhg1_shuff_bar_res{b}.^2,2)) ;
    inhg2_shuff_bar_var(b) = mean(mean(inhg2_shuff_bar_res{b}.^2,2)) ;
end


% g correlations

for b=1:numBars ; % preallocate variables to improve memory usage and speed
    c1_plusConverge{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    c2_plusConverge{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    c1_minusConverge{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    c2_minusConverge{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    cExc_plusConverg_plusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    cInh_plusConverg_plusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;  
    cExc_minusConverg_plusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    cInh_minusConverg_plusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;   
    cExc_plusConverg_minusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    cInh_plusConverg_minusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    cExc_minusConverg_minusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    cInh_minusConverg_minusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ; 
    cExcInh_plusConverg_plusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    cInhExc_plusConverg_plusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    cExcInh_minusConverg_plusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    cInhExc_minusConverg_plusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    cExcInh_plusConverg_minusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    cInhExc_plusConverg_minusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    cExcInh_minusConverg_minusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    cInhExc_minusConverg_minusPairwise{b} = nan(numTrials-2*residualOption,barPnts*2-1) ; 
    excg1_sim_bar_res_ac{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    inhg1_sim_bar_res_ac{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    excg2_sim_bar_res_ac{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    inhg2_sim_bar_res_ac{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    excg1_shuff_bar_res_ac{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    inhg1_shuff_bar_res_ac{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    excg2_shuff_bar_res_ac{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    inhg2_shuff_bar_res_ac{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    excg2_sim_bar_res_circshift_ac{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    inhg2_sim_bar_res_circshift_ac{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    excg2_shuff_bar_res_circshift_ac{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
    inhg2_shuff_bar_res_circshift_ac{b} = nan(numTrials-2*residualOption,barPnts*2-1) ;
end
c1_plusConverge_peak = nan(numBars,numTrials-2*residualOption) ; 
c2_plusConverge_peak = nan(numBars,numTrials-2*residualOption) ;  
c1_minusConverge_peak = nan(numBars,numTrials-2*residualOption) ; 
c2_minusConverge_peak = nan(numBars,numTrials-2*residualOption) ; 
cExc_plusConverg_plusPairwise_peak = nan(numBars,numTrials-2*residualOption) ; 
cInh_plusConverg_plusPairwise_peak = nan(numBars,numTrials-2*residualOption) ;         
cExc_minusConverg_plusPairwise_peak = nan(numBars,numTrials-2*residualOption) ; 
cInh_minusConverg_plusPairwise_peak = nan(numBars,numTrials-2*residualOption) ;        
cExc_plusConverg_minusPairwise_peak = nan(numBars,numTrials-2*residualOption) ; 
cInh_plusConverg_minusPairwise_peak = nan(numBars,numTrials-2*residualOption) ;        
cExc_minusConverg_minusPairwise_peak = nan(numBars,numTrials-2*residualOption) ; 
cInh_minusConverg_minusPairwise_peak = nan(numBars,numTrials-2*residualOption) ;         
cExcInh_plusConverg_plusPairwise_peak = nan(numBars,numTrials-2*residualOption) ; 
cInhExc_plusConverg_plusPairwise_peak = nan(numBars,numTrials-2*residualOption) ;  
cExcInh_minusConverg_plusPairwise_peak = nan(numBars,numTrials-2*residualOption) ; 
cInhExc_minusConverg_plusPairwise_peak = nan(numBars,numTrials-2*residualOption) ;       
cExcInh_plusConverg_minusPairwise_peak = nan(numBars,numTrials-2*residualOption) ; 
cInhExc_plusConverg_minusPairwise_peak = nan(numBars,numTrials-2*residualOption) ; 
cExcInh_minusConverg_minusPairwise_peak = nan(numBars,numTrials-2*residualOption) ; 
cInhExc_minusConverg_minusPairwise_peak = nan(numBars,numTrials-2*residualOption) ;  
cSpikePred_plusConverg_plusPairwise = nan(numBars,numTrials-2*residualOption) ; 
cSpikePred_minusConverg_plusPairwise = nan(numBars,numTrials-2*residualOption) ; 
cSpikePred_plusConverg_minusPairwise = nan(numBars,numTrials-2*residualOption) ; 
cSpikePred_minusConverg_minusPairwise = nan(numBars,numTrials-2*residualOption) ;         
c1_plusConverge_mean = nan(numBars,barPnts*2-1) ; 
c2_plusConverge_mean = nan(numBars,barPnts*2-1) ;    
c1_minusConverge_mean = nan(numBars,barPnts*2-1) ; 
c2_minusConverge_mean = nan(numBars,barPnts*2-1) ;     
cExc_plusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cInh_plusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cExc_minusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cInh_minusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cExc_plusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cInh_plusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cExc_minusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cInh_minusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ;     
cExcInh_plusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cInhExc_plusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cExcInh_minusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cInhExc_minusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cExcInh_plusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cInhExc_plusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cExcInh_minusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ; 
cInhExc_minusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ;  


for b=1:numBars ;
    excg2_sim_bar_res_circshift{b} = circshift(excg2_sim_bar_res{b},[minusPairwiseShift,0]) ;
    inhg2_sim_bar_res_circshift{b} = circshift(inhg2_sim_bar_res{b},[minusPairwiseShift,0]) ;
    excg2_shuff_bar_res_circshift{b} = circshift(excg2_shuff_bar_res{b},[minusPairwiseShift,0]) ;
    inhg2_shuff_bar_res_circshift{b} = circshift(inhg2_shuff_bar_res{b},[minusPairwiseShift,0]) ;
    
    for a = 1:numTrials-2*residualOption ; % for each residual     
        
        % converging correlations unormalized
        c1_plusConverge{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),inhg1_sim_bar_res{b}(a,:),'biased') ;
        c2_plusConverge{b}(a,:) = xcorr(excg2_sim_bar_res{b}(a,:),inhg2_sim_bar_res{b}(a,:),'biased') ;
        
        c1_minusConverge{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),inhg1_shuff_bar_res{b}(a,:),'biased') ;
        c2_minusConverge{b}(a,:) = xcorr(excg2_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res{b}(a,:),'biased') ;
    

        % pairwise ee,ii correlations unormalized
        cExc_plusConverg_plusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),excg2_sim_bar_res{b}(a,:),'biased') ; % + converging + pairwise
        cInh_plusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res{b}(a,:),'biased') ;
        
        cExc_minusConverg_plusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res{b}(a,:),'biased') ; % -converging + pairwise
        cInh_minusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res{b}(a,:),'biased') ;
        
        cExc_plusConverg_minusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),excg2_sim_bar_res_circshift{b}(a,:),'biased') ; % + converging - pairwise
        cInh_plusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res_circshift{b}(a,:),'biased') ;
        
        cExc_minusConverg_minusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res_circshift{b}(a,:),'unbiased') ; % -converging - pairwise
        cInh_minusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res_circshift{b}(a,:),'unbiased') ; 
  
        % pairwise ei,ie correlations unormalized
        cExcInh_plusConverg_plusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res{b}(a,:),'biased') ; % + converging + pairwise
        cInhExc_plusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),excg2_sim_bar_res{b}(a,:),'biased') ;
        
        cExcInh_minusConverg_plusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res{b}(a,:),'biased') ; % -converging + pairwise
        cInhExc_minusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res{b}(a,:),'biased') ;
        
        cExcInh_plusConverg_minusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res_circshift{b}(a,:),'biased') ; % + converging - pairwise
        cInhExc_plusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),excg2_sim_bar_res_circshift{b}(a,:),'biased') ;
        
        cExcInh_minusConverg_minusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res_circshift{b}(a,:),'biased') ; % -converging - pairwise
        cInhExc_minusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res_circshift{b}(a,:),'biased') ;    
        
        % c peaks 
        c1_plusConverge_peak(b,a) = CCpeakFinder(c1_plusConverge{b}(a,:)) ;
        c2_plusConverge_peak(b,a) = CCpeakFinder(c2_plusConverge{b}(a,:)) ;
        
        c1_minusConverge_peak(b,a) = CCpeakFinder(c1_minusConverge{b}(a,:)) ;
        c2_minusConverge_peak(b,a) = CCpeakFinder(c2_minusConverge{b}(a,:)) ;
    
        cExc_plusConverg_plusPairwise_peak(b,a) = CCpeakFinder(cExc_plusConverg_plusPairwise{b}(a,:)) ; % + converging + pairwise
        cInh_plusConverg_plusPairwise_peak(b,a) = CCpeakFinder(cInh_plusConverg_plusPairwise{b}(a,:)) ;
        
        cExc_minusConverg_plusPairwise_peak(b,a) = CCpeakFinder(cExc_minusConverg_plusPairwise{b}(a,:)) ; % -converging + pairwise
        cInh_minusConverg_plusPairwise_peak(b,a) = CCpeakFinder(cInh_minusConverg_plusPairwise{b}(a,:)) ;
        
        cExc_plusConverg_minusPairwise_peak(b,a) = CCpeakFinder(cExc_plusConverg_minusPairwise{b}(a,:)) ; % + converging - pairwise
        cInh_plusConverg_minusPairwise_peak(b,a) = CCpeakFinder(cInh_plusConverg_minusPairwise{b}(a,:)) ;
        
        cExc_minusConverg_minusPairwise_peak(b,a) = CCpeakFinder(cExc_minusConverg_minusPairwise{b}(a,:)) ; % -converging - pairwise
        cInh_minusConverg_minusPairwise_peak(b,a) = CCpeakFinder(cInh_minusConverg_minusPairwise{b}(a,:)) ;   
        
        cExcInh_plusConverg_plusPairwise_peak(b,a) = CCpeakFinder(cExcInh_plusConverg_plusPairwise{b}(a,:)) ; % + converging + pairwise
        cInhExc_plusConverg_plusPairwise_peak(b,a) = CCpeakFinder(cInhExc_plusConverg_plusPairwise{b}(a,:)) ;
        
        cExcInh_minusConverg_plusPairwise_peak(b,a) = CCpeakFinder(cExcInh_minusConverg_plusPairwise{b}(a,:)) ; % -converging + pairwise
        cInhExc_minusConverg_plusPairwise_peak(b,a) = CCpeakFinder(cInhExc_minusConverg_plusPairwise{b}(a,:)) ;
        
        cExcInh_plusConverg_minusPairwise_peak(b,a) = CCpeakFinder(cExcInh_plusConverg_minusPairwise{b}(a,:)) ; % + converging - pairwise
        cInhExc_plusConverg_minusPairwise_peak(b,a) = CCpeakFinder(cInhExc_plusConverg_minusPairwise{b}(a,:)) ;
        
        cExcInh_minusConverg_minusPairwise_peak(b,a) = CCpeakFinder(cExcInh_minusConverg_minusPairwise{b}(a,:)) ; % -converging - pairwise
        cInhExc_minusConverg_minusPairwise_peak(b,a) = CCpeakFinder(cInhExc_minusConverg_minusPairwise{b}(a,:)) ;   
        
        % autocorrelations
        excg1_sim_bar_res_ac{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),'biased') ;
        inhg1_sim_bar_res_ac{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),'biased') ;
        excg2_sim_bar_res_ac{b}(a,:) = xcorr(excg2_sim_bar_res{b}(a,:),'biased') ;
        inhg2_sim_bar_res_ac{b}(a,:) = xcorr(inhg2_sim_bar_res{b}(a,:),'biased') ;
        excg1_shuff_bar_res_ac{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),'biased') ;
        inhg1_shuff_bar_res_ac{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),'biased') ;
        excg2_shuff_bar_res_ac{b}(a,:) = xcorr(excg2_shuff_bar_res{b}(a,:),'biased') ;
        inhg2_shuff_bar_res_ac{b}(a,:) = xcorr(inhg2_shuff_bar_res{b}(a,:),'biased') ;
        excg2_sim_bar_res_circshift_ac{b}(a,:) = xcorr(excg2_sim_bar_res_circshift{b}(a,:),'biased') ;
        inhg2_sim_bar_res_circshift_ac{b}(a,:) = xcorr(inhg2_sim_bar_res_circshift{b}(a,:),'biased') ;
        excg2_shuff_bar_res_circshift_ac{b}(a,:) = xcorr(excg2_shuff_bar_res_circshift{b}(a,:),'biased') ;
        inhg2_shuff_bar_res_circshift_ac{b}(a,:) = xcorr(inhg2_shuff_bar_res_circshift{b}(a,:),'biased') ;
        
 
    end

    % mean corrs
    c1_plusConverge_mean(b,:) = mean(c1_plusConverge{b}) ;
    c2_plusConverge_mean(b,:) = mean(c2_plusConverge{b}) ;
    
    c1_minusConverge_mean(b,:) = mean(c1_minusConverge{b}) ;
    c2_minusConverge_mean(b,:) = mean(c2_minusConverge{b}) ; 
    
    cExc_plusConverg_plusPairwise_mean(b,:) = mean(cExc_plusConverg_plusPairwise{b}) ;
    cInh_plusConverg_plusPairwise_mean(b,:) = mean(cInh_plusConverg_plusPairwise{b}) ;

    cExc_minusConverg_plusPairwise_mean(b,:) = mean(cExc_minusConverg_plusPairwise{b}) ;
    cInh_minusConverg_plusPairwise_mean(b,:) = mean(cInh_minusConverg_plusPairwise{b}) ;
    
    cExc_plusConverg_minusPairwise_mean(b,:) = mean(cExc_plusConverg_minusPairwise{b}) ; 
    cInh_plusConverg_minusPairwise_mean(b,:) = mean(cInh_plusConverg_minusPairwise{b}) ;

    cExc_minusConverg_minusPairwise_mean(b,:) = mean(cExc_minusConverg_minusPairwise{b}) ; 
    cInh_minusConverg_minusPairwise_mean(b,:) = mean(cInh_minusConverg_minusPairwise{b}) ;   
    
    cExcInh_plusConverg_plusPairwise_mean(b,:) = mean(cExcInh_plusConverg_plusPairwise{b}) ;
    cInhExc_plusConverg_plusPairwise_mean(b,:) = mean(cInhExc_plusConverg_plusPairwise{b}) ;

    cExcInh_minusConverg_plusPairwise_mean(b,:) = mean(cExcInh_minusConverg_plusPairwise{b}) ;
    cInhExc_minusConverg_plusPairwise_mean(b,:) = mean(cInhExc_minusConverg_plusPairwise{b}) ;

    cExcInh_plusConverg_minusPairwise_mean(b,:) = mean(cExcInh_plusConverg_minusPairwise{b}) ;
    cInhExc_plusConverg_minusPairwise_mean(b,:) = mean(cInhExc_plusConverg_minusPairwise{b}) ;

    cExcInh_minusConverg_minusPairwise_mean(b,:) = mean(cExcInh_minusConverg_minusPairwise{b}) ;
    cInhExc_minusConverg_minusPairwise_mean(b,:) = mean(cInhExc_minusConverg_minusPairwise{b}) ;  
    
    % peak of mean corrs
    c1_plusConverge_mean_peak(b) = CCpeakFinder(c1_plusConverge_mean(b,:)) ;
    c2_plusConverge_mean_peak(b) = CCpeakFinder(c2_plusConverge_mean(b,:)) ;
    
    c1_minusConverge_mean_peak(b) = CCpeakFinder(c1_minusConverge_mean(b,:)) ;
    c2_minusConverge_mean_peak(b) = CCpeakFinder(c2_minusConverge_mean(b,:)) ; 
    
    cExc_plusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(cExc_plusConverg_plusPairwise_mean(b,:)) ;
    cInh_plusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(cInh_plusConverg_plusPairwise_mean(b,:)) ;

    cExc_minusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(cExc_minusConverg_plusPairwise_mean(b,:)) ;
    cInh_minusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(cInh_minusConverg_plusPairwise_mean(b,:)) ;
    
    cExc_plusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(cExc_plusConverg_minusPairwise_mean(b,:)) ; 
    cInh_plusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(cInh_plusConverg_minusPairwise_mean(b,:)) ;

    cExc_minusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(cExc_minusConverg_minusPairwise_mean(b,:)) ; 
    cInh_minusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(cInh_minusConverg_minusPairwise_mean(b,:)) ;  
    
    cExcInh_plusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(cExcInh_plusConverg_plusPairwise_mean(b,:)) ;
    cInhExc_plusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(cInhExc_plusConverg_plusPairwise_mean(b,:)) ;

    cExcInh_minusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(cExcInh_minusConverg_plusPairwise_mean(b,:)) ;
    cInhExc_minusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(cInhExc_minusConverg_plusPairwise_mean(b,:)) ;

    cExcInh_plusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(cExcInh_plusConverg_minusPairwise_mean(b,:)) ;
    cInhExc_plusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(cInhExc_plusConverg_minusPairwise_mean(b,:)) ;

    cExcInh_minusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(cExcInh_minusConverg_minusPairwise_mean(b,:)) ;
    cInhExc_minusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(cInhExc_minusConverg_minusPairwise_mean(b,:)) ;  
    
    % autocorrelations means
    excg1_sim_bar_res_ac_mean(b,:) = mean(excg1_sim_bar_res_ac{b}) ;
    inhg1_sim_bar_res_ac_mean(b,:) = mean(inhg1_sim_bar_res_ac{b}) ;
    excg2_sim_bar_res_ac_mean(b,:) = mean(excg2_sim_bar_res_ac{b}) ;
    inhg2_sim_bar_res_ac_mean(b,:) = mean(inhg2_sim_bar_res_ac{b}) ;
    excg1_shuff_bar_res_ac_mean(b,:) = mean(excg1_shuff_bar_res_ac{b}) ;
    inhg1_shuff_bar_res_ac_mean(b,:) = mean(inhg1_shuff_bar_res_ac{b}) ;
    excg2_shuff_bar_res_ac_mean(b,:) = mean(excg2_shuff_bar_res_ac{b}) ;
    inhg2_shuff_bar_res_ac_mean(b,:) = mean(inhg2_shuff_bar_res_ac{b}) ;
    excg2_sim_bar_res_circshift_ac_mean(b,:) = mean(excg2_sim_bar_res_circshift_ac{b}) ;
    inhg2_sim_bar_res_circshift_ac_mean(b,:) = mean(inhg2_sim_bar_res_circshift_ac{b}) ;
    excg2_shuff_bar_res_circshift_ac_mean(b,:) = mean(excg2_shuff_bar_res_circshift_ac{b}) ;
    inhg2_shuff_bar_res_circshift_ac_mean(b,:) = mean(inhg2_shuff_bar_res_circshift_ac{b}) ;

    % autocorrelations mean peaks
    
    excg1_sim_bar_res_ac_mean_peak(b) = max(excg1_sim_bar_res_ac_mean(b,:)) ;
    inhg1_sim_bar_res_ac_mean_peak(b) = max(inhg1_sim_bar_res_ac_mean(b,:)) ;
    excg2_sim_bar_res_ac_mean_peak(b) = max(excg2_sim_bar_res_ac_mean(b,:)) ;
    inhg2_sim_bar_res_ac_mean_peak(b) = max(inhg2_sim_bar_res_ac_mean(b,:)) ;
    excg1_shuff_bar_res_ac_mean_peak(b) = max(excg1_shuff_bar_res_ac_mean(b,:)) ;
    inhg1_shuff_bar_res_ac_mean_peak(b) = max(inhg1_shuff_bar_res_ac_mean(b,:)) ;
    excg2_shuff_bar_res_ac_mean_peak(b) = max(excg2_shuff_bar_res_ac_mean(b,:)) ;
    inhg2_shuff_bar_res_ac_mean_peak(b) = max(inhg2_shuff_bar_res_ac_mean(b,:)) ;
    excg2_sim_bar_res_circshift_ac_mean_peak(b) = max(excg2_sim_bar_res_circshift_ac_mean(b,:)) ;
    inhg2_sim_bar_res_circshift_ac_mean_peak(b) = max(inhg2_sim_bar_res_circshift_ac_mean(b,:)) ;
    excg2_shuff_bar_res_circshift_ac_mean_peak(b) = max(excg2_shuff_bar_res_circshift_ac_mean(b,:)) ;
    inhg2_shuff_bar_res_circshift_ac_mean_peak(b) = max(inhg2_shuff_bar_res_circshift_ac_mean(b,:)) ;
    
end

% linear analytic prediction
IntegratedCCOption = 1 ;

if IntegratedCCOption == 1 ;
     
    for c = 1:length(cExc_plusConverg_plusPairwise_mean) ; % for every lag

        % analytical solution  
        % Pss = (Ve1e2 +a^2*Vi1i2 -a*Ve1i2 -a*Ve2i1)/sqrt((Ve1^2 +a^2Vi1^2 -2aVe1i1)*(Ve2^2 +a^2Vi2^2 -2aVe2i2))
        % a (alpha) is ratio of inh/exc
        Ve1e2 = [cExc_plusConverg_plusPairwise_mean(:,c),cExc_minusConverg_plusPairwise_mean(:,c),cExc_plusConverg_minusPairwise_mean(:,c),cExc_minusConverg_minusPairwise_mean(:,c)] ;
        Vi1i2 = [cInh_plusConverg_plusPairwise_mean(:,c),cInh_minusConverg_plusPairwise_mean(:,c),cInh_plusConverg_minusPairwise_mean(:,c),cInh_minusConverg_minusPairwise_mean(:,c)] ;
        Ve1i2 = [cExcInh_plusConverg_plusPairwise_mean(:,c),cExcInh_minusConverg_plusPairwise_mean(:,c),cExcInh_plusConverg_minusPairwise_mean(:,c),cExcInh_minusConverg_minusPairwise_mean(:,c)] ;
        Vi1e2 = [cInhExc_plusConverg_plusPairwise_mean(:,c),cInhExc_minusConverg_plusPairwise_mean(:,c),cInhExc_plusConverg_minusPairwise_mean(:,c),cInhExc_minusConverg_minusPairwise_mean(:,c)] ;

        Ve1i1 = [c1_plusConverge_mean(:,c),c1_minusConverge_mean(:,c),c1_plusConverge_mean(:,c),c1_minusConverge_mean(:,c)] ;
        Ve2i2 = [c2_plusConverge_mean(:,c),c2_minusConverge_mean(:,c),c2_plusConverge_mean(:,c),c2_minusConverge_mean(:,c)] ;
        
        Ve1 = [excg1_sim_bar_res_ac_mean(:,c),excg1_shuff_bar_res_ac_mean(:,c),excg1_sim_bar_res_ac_mean(:,c),excg1_shuff_bar_res_ac_mean(:,c)] ;
        Ve2 = [excg2_sim_bar_res_ac_mean(:,c),excg2_shuff_bar_res_ac_mean(:,c),excg2_sim_bar_res_ac_mean(:,c),excg2_shuff_bar_res_ac_mean(:,c)] ;
        Vi1 = [inhg1_sim_bar_res_ac_mean(:,c),inhg1_shuff_bar_res_ac_mean(:,c),inhg1_sim_bar_res_ac_mean(:,c),inhg1_shuff_bar_res_ac_mean(:,c)] ;
        Vi2 = [inhg2_sim_bar_res_ac_mean(:,c),inhg2_shuff_bar_res_ac_mean(:,c),inhg2_sim_bar_res_ac_mean(:,c),inhg2_shuff_bar_res_ac_mean(:,c)] ;


        numConditions = 4 ; 
        alpha =  abs((mean(Vthresh_opt)-InhRev)/(mean(Vthresh_opt)-0)) ; % alpha from LIF spike threshold
        for a=1:numBars; % for every bar
            for b=1:numConditions ; % for every set pairings   
                sn_covariance_predicted_cc{b}(a,c) = Ve1e2(a,b) + alpha^2*Vi1i2(a,b) - alpha*Ve1i2(a,b) - alpha*Vi1e2(a,b) ;
                sn_variance1_predicted_ac{b}(a,c) = Ve1(a,b) + alpha^2*Vi1(a,b) - 2*alpha*Ve1i1(a,b) ;
                sn_variance2_predicted_ac{b}(a,c) = Ve2(a,b) + alpha^2*Vi2(a,b) - 2*alpha*Ve2i2(a,b) ;
            end            
        end
    end

    for b=1:numConditions ; % for every set pairings sum accross lag with scale factor
        sn_covariance_predicted(:,b) = ((length(sn_covariance_predicted_cc{b})-1)/2)*sum(sn_covariance_predicted_cc{b},2) ;
        sn_variance1_predicted(:,b) = ((length(sn_variance1_predicted_ac{b})-1)/2)*sum(sn_variance1_predicted_ac{b},2) ;
        sn_variance2_predicted(:,b) = ((length(sn_variance2_predicted_ac{b})-1)/2)*sum(sn_variance2_predicted_ac{b},2) ;
    
    end
    sn_corrcoef_predicted = sn_covariance_predicted./sqrt(sn_variance1_predicted.*sn_variance2_predicted) ; 
    
    
else

    % analytical solution  
    % Pss = (Ve1e2 +a^2*Vi1i2 -a*Ve1i2 -a*Ve2i1)/sqrt((Ve1^2 +a^2Vi1^2 -2aVe1i1)*(Ve2^2 +a^2Vi2^2 -2aVe2i2))
    % a (alpha) is ratio of inh/exc
    Ve1e2 = [cExc_plusConverg_plusPairwise_mean_peak(:),cExc_minusConverg_plusPairwise_mean_peak(:),cExc_plusConverg_minusPairwise_mean_peak(:),cExc_minusConverg_minusPairwise_mean_peak(:)] ;
    Vi1i2 = [cInh_plusConverg_plusPairwise_mean_peak(:),cInh_minusConverg_plusPairwise_mean_peak(:),cInh_plusConverg_minusPairwise_mean_peak(:),cInh_minusConverg_minusPairwise_mean_peak(:)] ;
    Ve1i2 = [cExcInh_plusConverg_plusPairwise_mean_peak(:),cExcInh_minusConverg_plusPairwise_mean_peak(:),cExcInh_plusConverg_minusPairwise_mean_peak(:),cExcInh_minusConverg_minusPairwise_mean_peak(:)] ;
    Vi1e2 = [cInhExc_plusConverg_plusPairwise_mean_peak(:),cInhExc_minusConverg_plusPairwise_mean_peak(:),cInhExc_plusConverg_minusPairwise_mean_peak(:),cInhExc_minusConverg_minusPairwise_mean_peak(:)] ;

    Ve1 = [excg1_sim_bar_res_ac_mean_peak(:),excg1_shuff_bar_res_ac_mean_peak(:),excg1_sim_bar_res_ac_mean_peak(:),excg1_shuff_bar_res_ac_mean_peak(:)] ;
    Ve2 = [excg2_sim_bar_res_ac_mean_peak(:),excg2_shuff_bar_res_ac_mean_peak(:),excg2_sim_bar_res_ac_mean_peak(:),excg2_shuff_bar_res_ac_mean_peak(:)] ;
    Vi1 = [inhg1_sim_bar_res_ac_mean_peak(:),inhg1_shuff_bar_res_ac_mean_peak(:),inhg1_sim_bar_res_ac_mean_peak(:),inhg1_shuff_bar_res_ac_mean_peak(:)] ;
    Vi2 = [inhg2_sim_bar_res_ac_mean_peak(:),inhg2_shuff_bar_res_ac_mean_peak(:),inhg2_sim_bar_res_ac_mean_peak(:),inhg2_shuff_bar_res_ac_mean_peak(:)] ;
    Ve1i1 = [c1_plusConverge_mean_peak(:),c1_minusConverge_mean_peak(:),c1_plusConverge_mean_peak(:),c1_minusConverge_mean_peak(:)] ;
    Ve2i2 = [c2_plusConverge_mean_peak(:),c2_minusConverge_mean_peak(:),c2_plusConverge_mean_peak(:),c2_minusConverge_mean_peak(:)] ;

    numConditions = 4 ; 
    alpha =  abs((mean(Vthresh_opt)-InhRev)/(mean(Vthresh_opt)-0)) ; % alpha from LIF spike threshold
    for a=1:numBars; % for every bar
        for b=1:numConditions ; % for every set pairings   
            sn_covariance_predicted(a,b) = Ve1e2(a,b) + alpha^2*Vi1i2(a,b) - alpha*Ve1i2(a,b) - alpha*Vi1e2(a,b) ;
            sn_variance1_predicted(a,b) = Ve1(a,b) + alpha^2*Vi1(a,b) - 2*alpha*Ve1i1(a,b) ;
            sn_variance2_predicted(a,b) = Ve2(a,b) + alpha^2*Vi2(a,b) - 2*alpha*Ve2i2(a,b) ;
            sn_corrcoef_predicted(a,b) = sn_covariance_predicted(a,b)/sqrt(sn_variance1_predicted(a,b)*sn_variance2_predicted(a,b)) ;
        end            
    end
    
end

% linear anlaytic predicted spike correlation coefficients by bar and condition
for a = 1:numBars ;
    identifier = ['LapPkCoefsPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_corrcoef_predicted(a,1);
    
    identifier = ['LapPkCoefsMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_corrcoef_predicted(a,2);

    identifier = ['LapPkCoefsPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_corrcoef_predicted(a,3);
    
    identifier = ['LapPkCoefsMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_corrcoef_predicted(a,4);
    
    
    identifier = ['LapPkCovsPcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_covariance_predicted(a,1);
    
    identifier = ['LapPkCovsMcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_covariance_predicted(a,2);

    identifier = ['LapPkCovsPcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_covariance_predicted(a,3);
    
    identifier = ['LapPkCovsMcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_covariance_predicted(a,4);
    
 
    identifier = ['LapPkVar1Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_variance1_predicted(a,1);
    
    identifier = ['LapPkVar1Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_variance1_predicted(a,2);
    
    identifier = ['LapPkVar2Pc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_variance2_predicted(a,1);
    
    identifier = ['LapPkVar2Mc',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_variance2_predicted(a,2);
end

% % input covariance for each bar and condition
% for a = 1:numBars ;
%     % +ei +p
%     identifier = ['CovPkE1E2PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1excg2_sim_bar_cov(a) ;
% 
%     identifier = ['CovPkI1I2PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1inhg2_sim_bar_cov(a) ;
%     
%     identifier = ['CovPkE1I2PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1inhg2_sim_bar_cov(a) ;    
% 
%     identifier = ['CovPkE2I1PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1excg2_sim_bar_cov(a) ;
%    
%     identifier = ['CovPkE1I1PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1inhg1_sim_bar_cov(a) ;
%     
%     identifier = ['CovPkE2I2PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg2inhg2_sim_bar_cov(a) ;
%     
%     % -ei +p
%     identifier = ['CovPkE1E2McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1excg2_shuff_bar_cov(a) ;
% 
%     identifier = ['CovPkI1I2McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1inhg2_shuff_bar_cov(a) ;
%     
%     identifier = ['CovPkE1I2McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1inhg2_shuff_bar_cov(a) ;    
% 
%     identifier = ['CovPkE2I1McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1excg2_shuff_bar_cov(a) ;
%    
%     identifier = ['CovPkE1I1McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1inhg1_shuff_bar_cov(a) ;
%     
%     identifier = ['CovPkE2I2McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg2inhg2_shuff_bar_cov(a) ;
%     
%     % -ei -p
%     identifier = ['CovPkE1E2McMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1excg2Shift_shuff_bar_cov(a) ;
% 
%     identifier = ['CovPkI1I2McMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1inhg2Shift_shuff_bar_cov(a) ;
%     
%     identifier = ['CovPkE1I2McMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1inhg2Shift_shuff_bar_cov(a) ;    
% 
%     identifier = ['CovPkE2I1McMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1excg2Shift_shuff_bar_cov(a) ;
%    
%     % +ei -p
%     identifier = ['CovPkE1E2PcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1excg2Shift_sim_bar_cov(a) ;
% 
%     identifier = ['CovPkI1I2PcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1inhg2Shift_sim_bar_cov(a) ;
%     
%     identifier = ['CovPkE1I2PcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = excg1inhg2Shift_sim_bar_cov(a) ;    
% 
%     identifier = ['CovPkE2I1PcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = inhg1excg2Shift_sim_bar_cov(a) ;
% 
% end
%}
% make room in memmory
clearvars -except ForIgor Parameters id A excg* inhg* sn_corrCoefs num* prePnts barPnts minusPairwiseShift SI time* residualOption Vthresh_opt InhRev


%% noise correlation coefs in presented conductances  
%{
% g correlations
%
for b=1:numBars ; % preallocate variables to improve memory usage and speed
    cc1_plusConverge{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    cc2_plusConverge{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    cc1_minusConverge{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    cc2_minusConverge{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccExc_plusConverg_plusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccInh_plusConverg_plusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccExc_minusConverg_plusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccInh_minusConverg_plusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccExc_plusConverg_minusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccInh_plusConverg_minusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccExc_minusConverg_minusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccInh_minusConverg_minusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;  
    ccExcInh_plusConverg_plusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccInhExc_plusConverg_plusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;     
    ccExcInh_minusConverg_plusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccInhExc_minusConverg_plusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccExcInh_plusConverg_minusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccInhExc_plusConverg_minusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccExcInh_minusConverg_minusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;
    ccInhExc_minusConverg_minusPairwise{b} = nan(numTrials-residualOption*2,barPnts*2-1) ;  
end
cc1_plusConverge_peak = nan(numBars,numTrials-residualOption*2) ;
cc2_plusConverge_peak = nan(numBars,numTrials-residualOption*2) ;        
cc1_minusConverge_peak = nan(numBars,numTrials-residualOption*2) ;
cc2_minusConverge_peak = nan(numBars,numTrials-residualOption*2) ;
ccExc_plusConverg_plusPairwise_peak = nan(numBars,numTrials-residualOption*2) ;
ccInh_plusConverg_plusPairwise_peak = nan(numBars,numTrials-residualOption*2) ;       
ccExc_minusConverg_plusPairwise_peak = nan(numBars,numTrials-residualOption*2) ;
ccInh_minusConverg_plusPairwise_peak = nan(numBars,numTrials-residualOption*2) ;        
ccExc_plusConverg_minusPairwise_peak = nan(numBars,numTrials-residualOption*2) ;
ccInh_plusConverg_minusPairwise_peak = nan(numBars,numTrials-residualOption*2) ;       
ccExc_minusConverg_minusPairwise_peak = nan(numBars,numTrials-residualOption*2) ;
ccInh_minusConverg_minusPairwise_peak = nan(numBars,numTrials-residualOption*2) ;      
ccExcInh_plusConverg_plusPairwise_peak = nan(numBars,numTrials-residualOption*2) ; 
ccInhExc_plusConverg_plusPairwise_peak = nan(numBars,numTrials-residualOption*2) ; 
ccExcInh_minusConverg_plusPairwise_peak = nan(numBars,numTrials-residualOption*2) ; 
ccInhExc_minusConverg_plusPairwise_peak = nan(numBars,numTrials-residualOption*2) ;         
ccExcInh_plusConverg_minusPairwise_peak = nan(numBars,numTrials-residualOption*2) ; 
ccInhExc_plusConverg_minusPairwise_peak = nan(numBars,numTrials-residualOption*2) ;        
ccExcInh_minusConverg_minusPairwise_peak = nan(numBars,numTrials-residualOption*2) ; 
ccInhExc_minusConverg_minusPairwise_peak = nan(numBars,numTrials-residualOption*2) ;      
cc1_plusConverge_mean = nan(numBars,barPnts*2-1) ; 
cc2_plusConverge_mean = nan(numBars,barPnts*2-1) ;    
cc1_minusConverge_mean = nan(numBars,barPnts*2-1) ; 
cc2_minusConverge_mean = nan(numBars,barPnts*2-1) ;    
ccExc_plusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
ccInh_plusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
ccExc_minusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
ccInh_minusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ;    
ccExc_plusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ; 
ccInh_plusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ; 
ccExc_minusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ; 
ccInh_minusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ;   
ccExcInh_plusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
ccInhExc_plusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
ccExcInh_minusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
ccInhExc_minusConverg_plusPairwise_mean = nan(numBars,barPnts*2-1) ; 
ccExcInh_plusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ; 
ccInhExc_plusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ; 
ccExcInh_minusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ; 
ccInhExc_minusConverg_minusPairwise_mean = nan(numBars,barPnts*2-1) ;  


for b=1:numBars ;
    excg2_sim_bar_res_circshift{b} = circshift(excg2_sim_bar_res{b},[minusPairwiseShift,0]) ;
    inhg2_sim_bar_res_circshift{b} = circshift(inhg2_sim_bar_res{b},[minusPairwiseShift,0]) ;
    excg2_shuff_bar_res_circshift{b} = circshift(excg2_shuff_bar_res{b},[minusPairwiseShift,0]) ;
    inhg2_shuff_bar_res_circshift{b} = circshift(inhg2_shuff_bar_res{b},[minusPairwiseShift,0]) ;
    
    for a = 1:numTrials-residualOption*2 ; % for each residual
        
        % converging correlations coef
        cc1_plusConverge{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),inhg1_sim_bar_res{b}(a,:),'coeff') ;
        cc2_plusConverge{b}(a,:) = xcorr(excg2_sim_bar_res{b}(a,:),inhg2_sim_bar_res{b}(a,:),'coeff') ;
        
        cc1_minusConverge{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),inhg1_shuff_bar_res{b}(a,:),'coeff') ;
        cc2_minusConverge{b}(a,:) = xcorr(excg2_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res{b}(a,:),'coeff') ;
        
    
        % pairwise ee,ii correlations coefs
        ccExc_plusConverg_plusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),excg2_sim_bar_res{b}(a,:),'coeff') ; % + converging + pairwise
        ccInh_plusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res{b}(a,:),'coeff') ;
        
        ccExc_minusConverg_plusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res{b}(a,:),'coeff') ; % -converging + pairwise
        ccInh_minusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res{b}(a,:),'coeff') ;
        
        ccExc_plusConverg_minusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),excg2_sim_bar_res_circshift{b}(a,:),'coeff') ; % + converging - pairwise
        ccInh_plusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res_circshift{b}(a,:),'coeff') ;
        
        ccExc_minusConverg_minusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res_circshift{b}(a,:),'coeff') ; % -converging - pairwise
        ccInh_minusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res_circshift{b}(a,:),'coeff') ;  

        % pairwise ei,ie correlations coefs
        ccExcInh_plusConverg_plusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res{b}(a,:),'coeff') ; % + converging + pairwise
        ccInhExc_plusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),excg2_sim_bar_res{b}(a,:),'coeff') ;
        
        ccExcInh_minusConverg_plusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res{b}(a,:),'coeff') ; % -converging + pairwise
        ccInhExc_minusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res{b}(a,:),'coeff') ;
        
        ccExcInh_plusConverg_minusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res_circshift{b}(a,:),'coeff') ; % + converging - pairwise
        ccInhExc_plusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),excg2_sim_bar_res_circshift{b}(a,:),'coeff') ;
        
        ccExcInh_minusConverg_minusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res_circshift{b}(a,:),'coeff') ; % -converging - pairwise
        ccInhExc_minusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res_circshift{b}(a,:),'coeff') ;  
  
        
        % cc peaks 
        cc1_plusConverge_peak(b,a) = CCpeakFinder(cc1_plusConverge{b}(a,:)) ;
        cc2_plusConverge_peak(b,a) = CCpeakFinder(cc2_plusConverge{b}(a,:)) ;
        
        cc1_minusConverge_peak(b,a) = CCpeakFinder(cc1_minusConverge{b}(a,:)) ;
        cc2_minusConverge_peak(b,a) = CCpeakFinder(cc2_minusConverge{b}(a,:)) ;
    
        ccExc_plusConverg_plusPairwise_peak(b,a) = CCpeakFinder(ccExc_plusConverg_plusPairwise{b}(a,:)) ; % + converging + pairwise
        ccInh_plusConverg_plusPairwise_peak(b,a) = CCpeakFinder(ccInh_plusConverg_plusPairwise{b}(a,:)) ;
        
        ccExc_minusConverg_plusPairwise_peak(b,a) = CCpeakFinder(ccExc_minusConverg_plusPairwise{b}(a,:)) ; % -converging + pairwise
        ccInh_minusConverg_plusPairwise_peak(b,a) = CCpeakFinder(ccInh_minusConverg_plusPairwise{b}(a,:)) ;
        
        ccExc_plusConverg_minusPairwise_peak(b,a) = CCpeakFinder(ccExc_plusConverg_minusPairwise{b}(a,:)) ; % + converging - pairwise
        ccInh_plusConverg_minusPairwise_peak(b,a) = CCpeakFinder(ccInh_plusConverg_minusPairwise{b}(a,:)) ;
        
        ccExc_minusConverg_minusPairwise_peak(b,a) = CCpeakFinder(ccExc_minusConverg_minusPairwise{b}(a,:)) ; % -converging - pairwise
        ccInh_minusConverg_minusPairwise_peak(b,a) = CCpeakFinder(ccInh_minusConverg_minusPairwise{b}(a,:)) ;   
        
        ccExcInh_plusConverg_plusPairwise_peak(b,a) = CCpeakFinder(ccExcInh_plusConverg_plusPairwise{b}(a,:)) ; % + converging + pairwise
        ccInhExc_plusConverg_plusPairwise_peak(b,a) = CCpeakFinder(ccInhExc_plusConverg_plusPairwise{b}(a,:)) ;
        
        ccExcInh_minusConverg_plusPairwise_peak(b,a) = CCpeakFinder(ccExcInh_minusConverg_plusPairwise{b}(a,:)) ; % -converging + pairwise
        ccInhExc_minusConverg_plusPairwise_peak(b,a) = CCpeakFinder(ccInhExc_minusConverg_plusPairwise{b}(a,:)) ;
        
        ccExcInh_plusConverg_minusPairwise_peak(b,a) = CCpeakFinder(ccExcInh_plusConverg_minusPairwise{b}(a,:)) ; % + converging - pairwise
        ccInhExc_plusConverg_minusPairwise_peak(b,a) = CCpeakFinder(ccInhExc_plusConverg_minusPairwise{b}(a,:)) ;
        
        ccExcInh_minusConverg_minusPairwise_peak(b,a) = CCpeakFinder(ccExcInh_minusConverg_minusPairwise{b}(a,:)) ; % -converging - pairwise
        ccInhExc_minusConverg_minusPairwise_peak(b,a) = CCpeakFinder(ccInhExc_minusConverg_minusPairwise{b}(a,:)) ;   
        
 
    end
    % mean corr coefs
    cc1_plusConverge_mean(b,:) = mean(cc1_plusConverge{b}) ;
    cc2_plusConverge_mean(b,:) = mean(cc2_plusConverge{b}) ;
    
    cc1_minusConverge_mean(b,:) = mean(cc1_minusConverge{b}) ;
    cc2_minusConverge_mean(b,:) = mean(cc2_minusConverge{b}) ; 
    
    ccExc_plusConverg_plusPairwise_mean(b,:) = mean(ccExc_plusConverg_plusPairwise{b}) ;
    ccInh_plusConverg_plusPairwise_mean(b,:) = mean(ccInh_plusConverg_plusPairwise{b}) ;

    ccExc_minusConverg_plusPairwise_mean(b,:) = mean(ccExc_minusConverg_plusPairwise{b}) ;
    ccInh_minusConverg_plusPairwise_mean(b,:) = mean(ccInh_minusConverg_plusPairwise{b}) ;
    
    ccExc_plusConverg_minusPairwise_mean(b,:) = mean(ccExc_plusConverg_minusPairwise{b}) ; 
    ccInh_plusConverg_minusPairwise_mean(b,:) = mean(ccInh_plusConverg_minusPairwise{b}) ;

    ccExc_minusConverg_minusPairwise_mean(b,:) = mean(ccExc_minusConverg_minusPairwise{b}) ; 
    ccInh_minusConverg_minusPairwise_mean(b,:) = mean(ccInh_minusConverg_minusPairwise{b}) ;   
    
    ccExcInh_plusConverg_plusPairwise_mean(b,:) = mean(ccExcInh_plusConverg_plusPairwise{b}) ;
    ccInhExc_plusConverg_plusPairwise_mean(b,:) = mean(ccInhExc_plusConverg_plusPairwise{b}) ;

    ccExcInh_minusConverg_plusPairwise_mean(b,:) = mean(ccExcInh_minusConverg_plusPairwise{b}) ;
    ccInhExc_minusConverg_plusPairwise_mean(b,:) = mean(ccInhExc_minusConverg_plusPairwise{b}) ;

    ccExcInh_plusConverg_minusPairwise_mean(b,:) = mean(ccExcInh_plusConverg_minusPairwise{b}) ;
    ccInhExc_plusConverg_minusPairwise_mean(b,:) = mean(ccInhExc_plusConverg_minusPairwise{b}) ;

    ccExcInh_minusConverg_minusPairwise_mean(b,:) = mean(ccExcInh_minusConverg_minusPairwise{b}) ;
    ccInhExc_minusConverg_minusPairwise_mean(b,:) = mean(ccInhExc_minusConverg_minusPairwise{b}) ;  

    
    % peak of mean corr coefs
    cc1_plusConverge_mean_peak(b) = CCpeakFinder(cc1_plusConverge_mean(b,:)) ;
    cc2_plusConverge_mean_peak(b) = CCpeakFinder(cc2_plusConverge_mean(b,:)) ;
    
    cc1_minusConverge_mean_peak(b) = CCpeakFinder(cc1_minusConverge_mean(b,:)) ;
    cc2_minusConverge_mean_peak(b) = CCpeakFinder(cc2_minusConverge_mean(b,:)) ; 
    
    ccExc_plusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(ccExc_plusConverg_plusPairwise_mean(b,:)) ;
    ccInh_plusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(ccInh_plusConverg_plusPairwise_mean(b,:)) ;

    ccExc_minusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(ccExc_minusConverg_plusPairwise_mean(b,:)) ;
    ccInh_minusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(ccInh_minusConverg_plusPairwise_mean(b,:)) ;
    
    ccExc_plusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(ccExc_plusConverg_minusPairwise_mean(b,:)) ; 
    ccInh_plusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(ccInh_plusConverg_minusPairwise_mean(b,:)) ;

    ccExc_minusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(ccExc_minusConverg_minusPairwise_mean(b,:)) ; 
    ccInh_minusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(ccInh_minusConverg_minusPairwise_mean(b,:)) ;  
    
    ccExcInh_plusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(ccExcInh_plusConverg_plusPairwise_mean(b,:)) ;
    ccInhExc_plusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(ccInhExc_plusConverg_plusPairwise_mean(b,:)) ;

    ccExcInh_minusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(ccExcInh_minusConverg_plusPairwise_mean(b,:)) ;
    ccInhExc_minusConverg_plusPairwise_mean_peak(b) = CCpeakFinder(ccInhExc_minusConverg_plusPairwise_mean(b,:)) ;

    ccExcInh_plusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(ccExcInh_plusConverg_minusPairwise_mean(b,:)) ;
    ccInhExc_plusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(ccInhExc_plusConverg_minusPairwise_mean(b,:)) ;

    ccExcInh_minusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(ccExcInh_minusConverg_minusPairwise_mean(b,:)) ;
    ccInhExc_minusConverg_minusPairwise_mean_peak(b) = CCpeakFinder(ccInhExc_minusConverg_minusPairwise_mean(b,:)) ;  
    
    
end


% figures
%{
figure % g converging correlation coefs
for a=1:numBars ;
    subplot(2,numBars,a)
    plot(time_cc,cc1_plusConverge_mean(a,:),'k')
    hold on
    plot(time_cc,cc1_minusConverge_mean(a,:),'g')
    plot(time_cc,ones(1,length(time_cc))*.3,'--')
    xlim([-1 1])
    ylim([-.5 1])
    
    if a==1;
        ylabel('cell 1')
    end
    
    
    subplot(2,numBars,a+numBars)
    plot(time_cc,cc2_plusConverge_mean(a,:),'k')
    hold on
    plot(time_cc,cc2_minusConverge_mean(a,:),'g')
    plot(time_cc,ones(1,length(time_cc))*.3,'--')
    xlim([-1 1])
    ylim([-.5 1])
    
    if a==1;
        ylabel('cell 2')
    end
end

pause
close

figure % g pairwise correlation coefs
for a=1:numBars ;
    subplot(2,numBars,a)
    plot(time_cc,ccInh_plusConverg_plusPairwise_mean(a,:),'r')
    hold on
    plot(time_cc,ccInh_minusConverg_plusPairwise_mean(a,:),'r--')
    plot(time_cc,ccInh_plusConverg_minusPairwise_mean(a,:),'g')
    plot(time_cc,ccInh_minusConverg_minusPairwise_mean(a,:),'g--') 
    plot(time_cc,ones(1,length(time_cc))*.3,'--')
    xlim([-1 1])
    ylim([-.5 1])
    
    if a==1;
        ylabel('Inh')
    end
    
    subplot(2,numBars,numBars+a)
    plot(time_cc,ccExc_plusConverg_plusPairwise_mean(a,:),'b')
    hold on
    plot(time_cc,ccExc_minusConverg_plusPairwise_mean(a,:),'b--')
    plot(time_cc,ccExc_plusConverg_minusPairwise_mean(a,:),'g')
    plot(time_cc,ccExc_minusConverg_minusPairwise_mean(a,:),'g--')
    plot(time_cc,ones(1,length(time_cc))*.3,'--')
    xlim([-1 1])
    ylim([-.5 1])
    
    if a==1;
        ylabel('Exc')
    end
    
end

pause
close


figure % g pairwise e,i correlations
for a=1:numBars ;
    subplot(2,numBars,a)
    plot(time_cc,ccExcInh_plusConverg_plusPairwise_mean(a,:),'r')
    hold on
    plot(time_cc,ccExcInh_minusConverg_plusPairwise_mean(a,:),'r--')
    plot(time_cc,ccExcInh_plusConverg_minusPairwise_mean(a,:),'g')
    plot(time_cc,ccExcInh_minusConverg_minusPairwise_mean(a,:),'g--') 
    xlim([-1 1])
    ylim([-.5 1])
    
    if a==1;
        ylabel('Exc1Inh2')
    end
    
    subplot(2,numBars,numBars+a)
    plot(time_cc,ccInhExc_plusConverg_plusPairwise_mean(a,:),'r')
    hold on
    plot(time_cc,ccInhExc_minusConverg_plusPairwise_mean(a,:),'r--')
    plot(time_cc,ccInhExc_plusConverg_minusPairwise_mean(a,:),'g')
    plot(time_cc,ccInhExc_minusConverg_minusPairwise_mean(a,:),'g--') 
    xlim([-1 1])
    ylim([-.5 1])
    
    if a==1;
        ylabel('Inh1Exc2')
    end
end

pause
close


figure % optimal linear wieghting of input correlations
subplot(1,2,1)
plot(sn_corrCoefs(:), Pss_estimate_sn(:),'o')
hold on
plot(sT_corrCoefs(:), Pss_estimate_st(:),'+')
xlabel('actual')
ylabel('predicted')

subplot(1,2,2)
plot(PrinEigVec_sn,'o')
hold on
plot(PrinEigVec_st,'+')
legend('sn','st')
title('EigenVectors')

set(gca,'XTick',[1:5],'XTickLabel',{'Pee','Pii','Pei','Pie','Pss'})

pause
close


figure % peaks of corr coefs in each condition
for a = 1:numBars ;
    subplot(1,numBars,a)
    plot([1:4],Pee(a,:),'b-')
    hold on
    plot([1:4],Pii(a,:),'r-')
    plot([1:4],Pei(a,:),'y-')
    plot([1:4],Pie(a,:),'c-')
    
    legend('Pee','Pii','Pei','Pie')
    set(gca,'XTick',[1:4],'XTickLabel',{'+c+p','-c+p','+c-p','-c-p'})
    xlabel('condition')
    ylabel('input corr coefs')
end

pause
close
%}

% FOR IGOR

% input correlations coefficients for each bar and condition
for a = 1:numBars ;
    % +ei +p
    identifier = ['CCpeakE1E2PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccExc_plusConverg_plusPairwise_mean_peak(a);

    identifier = ['CCpeakI1I2PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccInh_plusConverg_plusPairwise_mean_peak(a);
    
    identifier = ['CCpeakE1I2PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccExcInh_plusConverg_plusPairwise_mean_peak(a);    

    identifier = ['CCpeakE2I1PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccInhExc_plusConverg_plusPairwise_mean_peak(a);
   
    identifier = ['CCpeakE1I1PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = cc1_plusConverge_mean_peak(a) ;
    
    identifier = ['CCpeakE2I2PcPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = cc2_plusConverge_mean_peak(a) ;
    
    % -ei +p
    identifier = ['CCpeakE1E2McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccExc_minusConverg_plusPairwise_mean_peak(a);

    identifier = ['CCpeakI1I2McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccInh_minusConverg_plusPairwise_mean_peak(a);
    
    identifier = ['CCpeakE1I2McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccExcInh_minusConverg_plusPairwise_mean_peak(a);    

    identifier = ['CCpeakE2I1McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccInhExc_minusConverg_plusPairwise_mean_peak(a);
   
    identifier = ['CCpeakE1I1McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = cc1_minusConverge_mean_peak(a) ;
    
    identifier = ['CCpeakE2I2McPp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = cc2_minusConverge_mean_peak(a) ;
    
    % -ei -p
    identifier = ['CCpeakE1E2McMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccExc_minusConverg_minusPairwise_mean_peak(a);

    identifier = ['CCpeakI1I2McMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccInh_minusConverg_minusPairwise_mean_peak(a);
    
    identifier = ['CCpeakE1I2McMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccExcInh_minusConverg_minusPairwise_mean_peak(a);    

    identifier = ['CCpeakE2I1McMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccInhExc_minusConverg_minusPairwise_mean_peak(a);
   
    % +ei -p
    identifier = ['CCpeakE1E2PcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccExc_plusConverg_minusPairwise_mean_peak(a);

    identifier = ['CCpeakI1I2PcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccInh_plusConverg_minusPairwise_mean_peak(a);
    
    identifier = ['CCpeakE1I2PcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccExcInh_plusConverg_minusPairwise_mean_peak(a);    

    identifier = ['CCpeakE2I1PcMp',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccInhExc_plusConverg_minusPairwise_mean_peak(a);
   

end


% % mean sim conductances
% identifier = ['GdcExc1Mean','id','cell',num2str(A)] ; % exc1
% ForIgor.(identifier) = excg1_sim_mean ; 
% 
% identifier = ['GdcInh1Mean','id','cell',num2str(A)] ; % inh1
% ForIgor.(identifier) = inhg1_sim_mean ; 
% 
% identifier = ['GdcExc2Mean','id','cell',num2str(A)] ; % exc2
% ForIgor.(identifier) = excg2_sim_mean ; 
% 
% identifier = ['GdcInh2Mean','id','cell',num2str(A)] ; % inh2
% ForIgor.(identifier) = inhg2_sim_mean ; 
% 
% % variance sim conductances
% identifier = ['GdcExc1Var','id','cell',num2str(A)] ; % exc1
% ForIgor.(identifier) = excg1_sim_var ; 
% 
% identifier = ['GdcInh1Var','id','cell',num2str(A)] ; % inh1
% ForIgor.(identifier) = inhg1_sim_var ; 
% 
% identifier = ['GdcExc2Var','id','cell',num2str(A)] ; % exc2
% ForIgor.(identifier) = excg2_sim_var ; 
% 
% identifier = ['GdcInh2Var','id','cell',num2str(A)] ; % inh2
% ForIgor.(identifier) = inhg2_sim_var ; 
% 
% % converging correlations
% identifier = ['Gdc_timecc','id','cell',num2str(A)] ; % exc1
% ForIgor.(identifier) = time_cc ; 
% 
% for a=1:numBars ;
% 
%     identifier = ['cc1_plusConvBar','id',num2str(a),'cell',num2str(A)] ; % cell 1 sim
%     ForIgor.(identifier) = cc1_plusConverge_mean(a,:) ; 
% 
%     identifier = ['cc2_plusConvBar','id',num2str(a),'cell',num2str(A)] ; % cell 2 sim
%     ForIgor.(identifier) = cc2_plusConverge_mean(a,:) ; 
%     
%     identifier = ['cc1_minusConvBar','id',num2str(a),'cell',num2str(A)] ; % cell 1 shuff
%     ForIgor.(identifier) = cc1_minusConverge_mean(a,:) ; 
% 
%     identifier = ['cc2_minusConvBar','id',num2str(a),'cell',num2str(A)] ; % cell 2 shuff
%     ForIgor.(identifier) = cc2_minusConverge_mean(a,:) ; 
% 
% end
%   


% % optimal linear weights and variances
% identifier = ['SNcorrWieghts','id',num2str(A)] ;
% ForIgor.(identifier) = PrinEigVec_sn' ;
% 
% identifier = ['STcorrWieghts','id',num2str(A)] ;
% ForIgor.(identifier) = PrinEigVec_st' ;
% 
% identifier = ['SNeigVariances','id',num2str(A)] ;
% ForIgor.(identifier) = sort(Variances_sn)' ;
% 
% identifier = ['STeigVariances','id',num2str(A)] ;
% ForIgor.(identifier) = sort(Variances_st)' ;

% alpha values


% identifier = ['snCorrsXlabel',num2str(a),'cell',num2str(A)] ;
% ForIgor.(identifier) = '+c+p','-c+p','+c-p','-c-p';   
    
%}   