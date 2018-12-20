function ForIgor = DCdsPairAnalyzer(Input,Parameters,id,A) ;

% this will analyze dynamic clamp data from on-off ds cells that were
% presented conductances recorded from pairs of ds cells during a moving
% bar.  Igor code presented conductances thusly in epochs 1:4
% 1) Cell1excTrial(t)+Cell1inhTrial(t) -->Sim 1
% 2) Cell2excTrial(t)+Cell2inhTrial(t) -->Sim 2
% 3) Cell1excTrial(t)+Cell1inhTrial(t+1) -->Shuffle 1
% 4) Cell2excTrial(t)+Cell2inhTrail(t+1) -->Shuffle 2
% looped over t (trial during condtance recording)
% THIS IS NOT IDEAL BECAUSE SIM AND SHUFFLE WILL NECESSARILY END WITH A
% COUPLE DIFFERENT INPUT CONDUCTANCES, but data is still useful.
%
% or thusly:
% 1) Cell1excTrial(t)+Cell1inhTrial(t) -->Sim 1
% 2) Cell2excTrial(t)+Cell2inhTrial(t) -->Sim 2
% 3) Cell1excTrial(t)+Cell1inhTrial(t+1) -->Shuffle 1
% 4) Cell2excTrial(t)+Cell2inhTrail(t+1) -->Shuffle 2
% 5) Cell1excTrial(t)+Cell1inhTrial(t) -->Sim 1
% 6) Cell2excTrial(t)+Cell2inhTrial(t) -->Sim 2
% 7) Cell1excTrial(t)+Cell1inhTrial(t-1) -->Shuffle 1
% 8) Cell2excTrial(t)+Cell2inhTrail(t-1) -->Shuffle 2
% THIS IS IDEAL

% JC 2/7/11

%% impact of changing pairwise e,i correlations
InhRev = -60 ; % reversal potential in igor code

epochs = str2num(Input(A).(id)) ; %sim and shuffled g

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);

% get data
[epochSize, error] = ITCGetEpochSize(epochs(1), fp) ; % preallocate variables to improve memory usage and speed 
data = nans(length(epochs),epochSize) ;
excg = nans(length(epochs),epochSize) ;
inhg = nans(length(epochs),epochSize) ;

for a = 1:length(epochs) ;
    [data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ; 
    [excg(a,:), inhg(a,:), error] = ITCReadEpochStmGClamp(epochs(a), 0, fp);
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
numBars = 3 ; % number of bars presented while recording conductances
prePnts = 2000 ; % pnts to avoid as g starts to be injected
barPnts = floor((length(data)-prePnts)/numBars) ; 

for a=1:numBars ; % preallocate variables to improve memory usage and speed
    spikeTrain1_sim_bar{a} = nans(numTrials,barPnts) ;
    spikeTrain2_sim_bar{a} = nans(numTrials,barPnts) ;
    spikeTrain1_shuff_bar{a} = nans(numTrials,barPnts) ;
    spikeTrain2_shuff_bar{a} = nans(numTrials,barPnts) ;
end

for a=1:numBars ;
    spikeTrain1_sim_bar{a} = spikeTrain1_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    spikeTrain2_sim_bar{a} = spikeTrain2_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    spikeTrain1_shuff_bar{a} = spikeTrain1_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    spikeTrain2_shuff_bar{a} = spikeTrain2_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
end

% spike number for each bar
sn1_sim = nans(numTrials,numBars) ; % preallocate variables to improve memory usage and speed
sn2_sim = nans(numTrials,numBars) ;
sn1_shuff = nans(numTrials,numBars) ;
sn2_shuff = nans(numTrials,numBars) ;
sn1_sim_mean = nans(1,numBars) ;
sn2_sim_mean = nans(1,numBars) ;
sn1_shuff_mean = nans(1,numBars) ;
sn2_shuff_mean = nans(1,numBars) ;

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

% smoothed spike trains
for a=1:numBars ; % preallocate variables to improve memory usage and speed
    spikeTrain1_smth_sim_bar{a} = nans(numTrials,barPnts) ;
    spikeTrain2_smth_sim_bar{a} = nans(numTrials,barPnts) ;
    spikeTrain1_smth_shuff_bar{a} = nans(numTrials,barPnts) ;
    spikeTrain2_smth_shuff_bar{a} = nans(numTrials,barPnts) ;    
end

smthFctr = .050 ; % sec
smthFctr = round(smthFctr/SI) ;
for a=1:numBars ;
    for b=1:numTrials ;
        spikeTrain1_smth_sim_bar{a}(b,:) = smooth(spikeTrain1_sim_bar{a}(b,:),smthFctr) ;
        spikeTrain2_smth_sim_bar{a}(b,:) = smooth(spikeTrain2_sim_bar{a}(b,:),smthFctr) ;
        spikeTrain1_smth_shuff_bar{a}(b,:) = smooth(spikeTrain1_shuff_bar{a}(b,:),smthFctr) ;
        spikeTrain2_smth_shuff_bar{a}(b,:) = smooth(spikeTrain2_shuff_bar{a}(b,:),smthFctr) ;    
    end
end

% spike number residuals (local correction)
sn1_sim_res = nans(numTrials-2,numBars) ; % preallocate variables to improve memory usage and speed
sn2_sim_res = nans(numTrials-2,numBars) ;
sn1_shuff_res = nans(numTrials-2,numBars) ;
sn2_shuff_res = nans(numTrials-2,numBars) ;

for a = 2:numTrials-1 ;
    sn1_sim_res(a-1,:) = sn1_sim(a,:) - (sn1_sim(a-1,:)+sn1_sim(a+1,:))/2 ;
    sn2_sim_res(a-1,:) = sn2_sim(a,:) - (sn2_sim(a-1,:)+sn2_sim(a+1,:))/2 ;
    sn1_shuff_res(a-1,:) = sn1_shuff(a,:) - (sn1_shuff(a-1,:)+sn1_shuff(a+1,:))/2 ;
    sn2_shuff_res(a-1,:) = sn2_shuff(a,:) - (sn2_shuff(a-1,:)+sn2_shuff(a+1,:))/2 ;
end

% smoothed spike train residuals
spikeTrain1_smth_sim_bar_res = cell(1,numBars) ; % preallocate variables to improve memory usage and speed
spikeTrain2_smth_sim_bar_res = cell(1,numBars) ;
spikeTrain1_smth_shuff_bar_res = cell(1,numBars) ;
spikeTrain2_smth_shuff_bar_res = cell(1,numBars) ;    

for a=1:numBars ;
    for b=2:numTrials-1  ;
        spikeTrain1_smth_sim_bar_res{a}(b-1,:) = spikeTrain1_smth_sim_bar{a}(b,:) - (spikeTrain1_smth_sim_bar{a}(b-1,:) + spikeTrain1_smth_sim_bar{a}(b+1,:))./2 ;
        spikeTrain2_smth_sim_bar_res{a}(b-1,:) = spikeTrain2_smth_sim_bar{a}(b,:) - (spikeTrain2_smth_sim_bar{a}(b-1,:) + spikeTrain2_smth_sim_bar{a}(b+1,:))./2 ;
        spikeTrain1_smth_shuff_bar_res{a}(b-1,:) = spikeTrain1_smth_shuff_bar{a}(b,:) - (spikeTrain1_smth_shuff_bar{a}(b-1,:) + spikeTrain1_smth_shuff_bar{a}(b+1,:))./2 ;
        spikeTrain2_smth_shuff_bar_res{a}(b-1,:) = spikeTrain2_smth_shuff_bar{a}(b,:) - (spikeTrain2_smth_shuff_bar{a}(b-1,:) + spikeTrain2_smth_shuff_bar{a}(b+1,:))./2 ;    
    end
end

% sn unbiased variance of spike number from residuals
sn1_sim_var = sum(sn1_sim_res.^2,1)./(size(sn1_sim_res,1)-1) ;
sn2_sim_var = sum(sn2_sim_res.^2,1)./(size(sn2_sim_res,1)-1) ;
sn1_shuff_var = sum(sn1_shuff_res.^2,1)./(size(sn1_shuff_res,1)-1) ;
sn2_shuff_var = sum(sn2_shuff_res.^2,1)./(size(sn2_shuff_res,1)-1) ;

% spike number correlation (linear) for each bar
lincoefsigMin = .05 ;
minusPairwiseShift = 3 ; % number of trials to shift to get rid of pairwise correlations
for a=1:numBars ;
    [Corr,p,lb,ub] =  corrcoef(sn1_sim_res(:,a),sn2_sim_res(:,a)) ; % + converging + pairwise correlations 
    sn_plusConverg_plusPairwise_lincoef(a) = Corr(1,2) ;
    sn_plusConverg_plusPairwise_lincoefsig(a) = p(1,2) ;
%     if sn_plusConverg_plusPairwise_lincoefsig(a)>lincoefsigMin ;
%         sn_plusConverg_plusPairwise_lincoef(a) = 0 ;
%     end
    sn_plusConverg_plusPairwise_std(a) = (Corr(1,2)-lb(1,2))*sqrt(length(sn1_sim_res(:,a)))/1.96 ;
    h=error_ellipse(Corr,[sn1_sim_mean(a),sn2_sim_mean(a)],'conf',.95) ; % confidence elipse
    sn_plusConverg_plusPairwise_ellipseX(a,:) = get(h,'XData') ;
    sn_plusConverg_plusPairwise_ellipseY(a,:) = get(h,'YData') ;
    
    [Corr,p,lb,ub] =  corrcoef(sn1_shuff_res(:,a),sn2_shuff_res(:,a)) ; % - converging + pairwise correlations
    sn_minusConverg_plusPairwise_lincoef(a) = Corr(1,2) ;
    sn_minusConverg_plusPairwise_lincoefsig(a) = p(1,2) ;
%     if sn_minusConverg_plusPairwise_lincoefsig(a)>lincoefsigMin ;
%         sn_minusConverg_plusPairwise_lincoef(a) = 0 ;
%     end
    sn_minusConverg_plusPairwise_std(a) = (Corr(1,2)-lb(1,2))*sqrt(length(sn1_sim_res(:,a)))/1.96 ;
    h=error_ellipse(Corr,[sn1_shuff_mean(a),sn2_shuff_mean(a)],'conf',.95) ; % confidence elipse
    sn_minusConverg_plusPairwise_ellipseX(a,:) = get(h,'XData') ;
    sn_minusConverg_plusPairwise_ellipseY(a,:) = get(h,'YData') ;
    
    [Corr,p,lb,ub] =  corrcoef(sn1_sim_res(:,a),circshift(sn2_sim_res(:,a),[minusPairwiseShift,0])) ; % + converging - pairwise correlations 
    sn_plusConverg_minusPairwise_lincoef(a) = Corr(1,2) ;
    sn_plusConverg_minusPairwise_lincoefsig(a) = p(1,2) ;
%     if sn_plusConverg_minusPairwise_lincoefsig(a)>lincoefsigMin ;
%         sn_plusConverg_minusPairwise_lincoef(a) = 0 ;
%     end
    sn_plusConverg_minusPairwise_std(a) = (Corr(1,2)-lb(1,2))*sqrt(length(sn1_sim_res(:,a)))/1.96 ;
    h=error_ellipse(Corr,[sn1_sim_mean(a),sn2_sim_mean(a)],'conf',.95) ; % confidence elipse
    sn_plusConverg_minusPairwise_ellipseX(a,:) = get(h,'XData') ;
    sn_plusConverg_minusPairwise_ellipseY(a,:) = get(h,'YData') ;
    
    [Corr,p,lb,ub] =  corrcoef(sn1_shuff_res(:,a),circshift(sn2_shuff_res(:,a),[minusPairwiseShift,0])) ; % - converging - pairwise correlations 
    sn_minusConverg_minusPairwise_lincoef(a) = Corr(1,2) ;
    sn_minusConverg_minusPairwise_lincoefsig(a) = p(1,2) ;
%     if sn_minusConverg_minusPairwise_lincoefsig(a)>lincoefsigMin ;
%         sn_minusConverg_minusPairwise_lincoef(a) = 0 ;
%     end
    sn_minusConverg_minusPairwise_std(a) = (Corr(1,2)-lb(1,2))*sqrt(length(sn1_sim_res(:,a)))/1.96 ;
    h=error_ellipse(Corr,[sn1_shuff_mean(a),sn2_shuff_mean(a)],'conf',.95) ; % confidence elipse
    sn_minusConverg_minusPairwise_ellipseX(a,:) = get(h,'XData') ;
    sn_minusConverg_minusPairwise_ellipseY(a,:) = get(h,'YData') ;
    
    sn_corrCoefs(a,:) = [sn_plusConverg_plusPairwise_lincoef(a),sn_minusConverg_plusPairwise_lincoef(a),...
    sn_plusConverg_minusPairwise_lincoef(a),sn_minusConverg_minusPairwise_lincoef(a)] ;
    
    sn_corrCoefsStd(a,:) = [sn_plusConverg_plusPairwise_std(a),sn_minusConverg_plusPairwise_std(a),...
    sn_plusConverg_minusPairwise_std(a),sn_minusConverg_minusPairwise_std(a)] ;
end
close,clear h

% smoothed spike train correlations
for a=1:numBars ; % preallocate variables to improve memory usage and speed
    spikeTrain2_smth_sim_bar_res_circshift{a} = nans(numTrials-2,barPnts) ; 
    spikeTrain2_smth_shuff_bar_res_circshift{a} = nans(numTrials-2,barPnts) ;        
    sT_plusConverg_plusPairwise_cc{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_minusConverg_plusPairwise_cc{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_plusConverg_minusPairwise_cc{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_minusConverg_minusPairwise_cc{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_plusConverg_plusPairwise_c{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_minusConverg_plusPairwise_c{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_plusConverg_minusPairwise_c{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_minusConverg_minusPairwise_c{a} = nans(numTrials-2,barPnts*2-1) ;
end
sT_plusConverg_plusPairwise_ccPeak = nans(numBars,numTrials-2) ;
sT_minusConverg_plusPairwise_ccPeak = nans(numBars,numTrials-2) ;
sT_plusConverg_minusPairwise_ccPeak = nans(numBars,numTrials-2) ;
sT_minusConverg_minusPairwise_ccPeak = nans(numBars,numTrials-2) ;
sT_plusConverg_plusPairwise_cPeak = nans(numBars,numTrials-2) ;
sT_minusConverg_plusPairwise_cPeak = nans(numBars,numTrials-2) ;
sT_plusConverg_minusPairwise_cPeak = nans(numBars,numTrials-2) ;
sT_minusConverg_minusPairwise_cPeak = nans(numBars,numTrials-2) ;
sT_plusConverg_plusPairwise_cc_mean = nans(numTrials-2,barPnts*2-1) ; 
sT_minusConverg_plusPairwise_cc_mean = nans(numTrials-2,barPnts*2-1) ;
sT_plusConverg_minusPairwise_cc_mean = nans(numTrials-2,barPnts*2-1) ;
sT_minusConverg_minusPairwise_cc_mean = nans(numTrials-2,barPnts*2-1) ;
sT_plusConverg_plusPairwise_c_mean = nans(numTrials-2,barPnts*2-1) ; 
sT_minusConverg_plusPairwise_c_mean = nans(numTrials-2,barPnts*2-1) ;
sT_plusConverg_minusPairwise_c_mean = nans(numTrials-2,barPnts*2-1) ;
sT_minusConverg_minusPairwise_c_mean = nans(numTrials-2,barPnts*2-1) ;
sT_plusConverg_plusPairwise_cc_mean_Peak = nans(1,numBars) ; 
sT_minusConverg_plusPairwise_cc_mean_Peak = nans(1,numBars) ;
sT_plusConverg_minusPairwise_cc_mean_Peak = nans(1,numBars) ;
sT_minusConverg_minusPairwise_cc_mean_Peak = nans(1,numBars) ;    
sT_plusConverg_plusPairwise_c_mean_Peak = nans(1,numBars) ; 
sT_minusConverg_plusPairwise_c_mean_Peak = nans(1,numBars) ;
sT_plusConverg_minusPairwise_c_mean_Peak = nans(1,numBars) ;
sT_minusConverg_minusPairwise_c_mean_Peak = nans(1,numBars) ; 


for a=1:numBars ;
    spikeTrain2_smth_sim_bar_res_circshift{a} = circshift(spikeTrain2_smth_sim_bar_res{a},[minusPairwiseShift,0]) ;
    spikeTrain2_smth_shuff_bar_res_circshift{a} = circshift(spikeTrain2_smth_shuff_bar_res{a},[minusPairwiseShift,0]) ;
    
    for b=1:numTrials-2  ;
        % corr coef
        sT_plusConverg_plusPairwise_cc{a}(b,:) = xcorr(spikeTrain1_smth_sim_bar_res{a}(b,:),spikeTrain2_smth_sim_bar_res{a}(b,:),'coef') ;
        sT_minusConverg_plusPairwise_cc{a}(b,:) = xcorr(spikeTrain1_smth_shuff_bar_res{a}(b,:),spikeTrain2_smth_shuff_bar_res{a}(b,:),'coef') ;
        sT_plusConverg_minusPairwise_cc{a}(b,:) = xcorr(spikeTrain1_smth_sim_bar_res{a}(b,:),spikeTrain2_smth_sim_bar_res_circshift{a}(b,:),'coef') ;
        sT_minusConverg_minusPairwise_cc{a}(b,:) = xcorr(spikeTrain1_smth_shuff_bar_res{a}(b,:),spikeTrain2_smth_shuff_bar_res_circshift{a}(b,:),'coef') ;
    
        % corr normalized by length only
        sT_plusConverg_plusPairwise_c{a}(b,:) = xcorr(spikeTrain1_smth_sim_bar_res{a}(b,:),spikeTrain2_smth_sim_bar_res{a}(b,:),'unbiased') ;
        sT_minusConverg_plusPairwise_c{a}(b,:) = xcorr(spikeTrain1_smth_shuff_bar_res{a}(b,:),spikeTrain2_smth_shuff_bar_res{a}(b,:),'unbiased') ;
        sT_plusConverg_minusPairwise_c{a}(b,:) = xcorr(spikeTrain1_smth_sim_bar_res{a}(b,:),spikeTrain2_smth_sim_bar_res_circshift{a}(b,:),'unbiased') ;
        sT_minusConverg_minusPairwise_c{a}(b,:) = xcorr(spikeTrain1_smth_shuff_bar_res{a}(b,:),spikeTrain2_smth_shuff_bar_res_circshift{a}(b,:),'unbiased') ;
        
        % peaks of cc
        sT_plusConverg_plusPairwise_ccPeak(a,b) = CCpeakFinder(sT_plusConverg_plusPairwise_cc{a}(b,:)) ;
        sT_minusConverg_plusPairwise_ccPeak(a,b) = CCpeakFinder(sT_minusConverg_plusPairwise_cc{a}(b,:)) ;
        sT_plusConverg_minusPairwise_ccPeak(a,b) = CCpeakFinder(sT_plusConverg_minusPairwise_cc{a}(b,:)) ;
        sT_minusConverg_minusPairwise_ccPeak(a,b) = CCpeakFinder(sT_minusConverg_minusPairwise_cc{a}(b,:)) ;
        
        % peaks of c
        sT_plusConverg_plusPairwise_cPeak(a,b) = CCpeakFinder(sT_plusConverg_plusPairwise_c{a}(b,:)) ;
        sT_minusConverg_plusPairwise_cPeak(a,b) = CCpeakFinder(sT_minusConverg_plusPairwise_c{a}(b,:)) ;
        sT_plusConverg_minusPairwise_cPeak(a,b) = CCpeakFinder(sT_plusConverg_minusPairwise_c{a}(b,:)) ;
        sT_minusConverg_minusPairwise_cPeak(a,b) = CCpeakFinder(sT_minusConverg_minusPairwise_c{a}(b,:)) ;
    end
    % mean corr coef
    sT_plusConverg_plusPairwise_cc_mean(a,:) = nanmean(sT_plusConverg_plusPairwise_cc{a}) ; % trials where residual = 0 have no cc (no residual = no noise = no ability to correlate noise)
    sT_minusConverg_plusPairwise_cc_mean(a,:) = nanmean(sT_minusConverg_plusPairwise_cc{a}) ;
    sT_plusConverg_minusPairwise_cc_mean(a,:) = nanmean(sT_plusConverg_minusPairwise_cc{a}) ;
    sT_minusConverg_minusPairwise_cc_mean(a,:) = nanmean(sT_minusConverg_minusPairwise_cc{a}) ;
    
    % mean corr 
    sT_plusConverg_plusPairwise_c_mean(a,:) = nanmean(sT_plusConverg_plusPairwise_c{a}) ; 
    sT_minusConverg_plusPairwise_c_mean(a,:) = nanmean(sT_minusConverg_plusPairwise_c{a}) ;
    sT_plusConverg_minusPairwise_c_mean(a,:) = nanmean(sT_plusConverg_minusPairwise_c{a}) ;
    sT_minusConverg_minusPairwise_c_mean(a,:) = nanmean(sT_minusConverg_minusPairwise_c{a}) ;
    
    % peak of mean corr coef
    sT_plusConverg_plusPairwise_cc_mean_Peak(a) = CCpeakFinder(sT_plusConverg_plusPairwise_cc_mean(a,:)) ; 
    sT_minusConverg_plusPairwise_cc_mean_Peak(a) = CCpeakFinder(sT_minusConverg_plusPairwise_cc_mean(a,:)) ;
    sT_plusConverg_minusPairwise_cc_mean_Peak(a) = CCpeakFinder(sT_plusConverg_minusPairwise_cc_mean(a,:)) ;
    sT_minusConverg_minusPairwise_cc_mean_Peak(a) = CCpeakFinder(sT_minusConverg_minusPairwise_cc_mean(a,:)) ;    
    
    % peak of mean corr
    sT_plusConverg_plusPairwise_c_mean_Peak(a) = CCpeakFinder(sT_plusConverg_plusPairwise_c_mean(a,:)) ; 
    sT_minusConverg_plusPairwise_c_mean_Peak(a) = CCpeakFinder(sT_minusConverg_plusPairwise_c_mean(a,:)) ;
    sT_plusConverg_minusPairwise_c_mean_Peak(a) = CCpeakFinder(sT_plusConverg_minusPairwise_c_mean(a,:)) ;
    sT_minusConverg_minusPairwise_c_mean_Peak(a) = CCpeakFinder(sT_minusConverg_minusPairwise_c_mean(a,:)) ;   
    
    % concatinated corr coefs
    sT_corrCoefs(a,:) = [sT_plusConverg_plusPairwise_cc_mean_Peak(a),sT_minusConverg_plusPairwise_cc_mean_Peak(a),...
    sT_plusConverg_minusPairwise_cc_mean_Peak(a),sT_minusConverg_minusPairwise_cc_mean_Peak(a)] ;

    % concatinated corrs
    sT_corr(a,:) = [sT_plusConverg_plusPairwise_c_mean_Peak(a),sT_minusConverg_plusPairwise_c_mean_Peak(a),...
    sT_plusConverg_minusPairwise_c_mean_Peak(a),sT_minusConverg_minusPairwise_c_mean_Peak(a)] ;
end

% conductances
excg1_sim = nans(numTrials,epochSize) ;  % preallocate variables to improve memory usage and speed
excg2_sim = nans(numTrials,epochSize) ;   
excg1_shuff = nans(numTrials,epochSize) ;  
excg2_shuff = nans(numTrials,epochSize) ;  
inhg1_sim = nans(numTrials,epochSize) ; 
inhg2_sim = nans(numTrials,epochSize) ;  
inhg1_shuff = nans(numTrials,epochSize) ;  
inhg2_shuff = nans(numTrials,epochSize) ;  

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
    excg1_sim_bar{a} = nans(numTrials,barPnts) ;  % preallocate variables to improve memory usage and speed 
    excg2_sim_bar{a} = nans(numTrials,barPnts) ;
    excg1_shuff_bar{a} = nans(numTrials,barPnts) ;
    excg2_shuff_bar{a} = nans(numTrials,barPnts) ;
    inhg1_sim_bar{a} = nans(numTrials,barPnts) ;
    inhg2_sim_bar{a} = nans(numTrials,barPnts) ;
    inhg1_shuff_bar{a} = nans(numTrials,barPnts) ;
    inhg2_shuff_bar{a} = nans(numTrials,barPnts) ;    
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
% figure % number of spikes per trial
% for a = 1:numBars ;
%     subplot(numBars,1,a)
%     plot(sn1_sim(:,a),'k')
%     hold on
%     plot(sn2_sim(:,a),'k:')
%     plot(sn1_shuff(:,a),'g')
%     plot(sn2_shuff(:,a),'g:')
% end
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
% figure % spike train noise correlations
% for a = 1:numBars ;
% subplot(1,numBars,a)
%     plot(time_cc,sT_plusConverg_plusPairwise_cc_mean(a,:),'k') 
%     hold on
%     plot(time_cc,sT_minusConverg_plusPairwise_cc_mean(a,:),'g')
%     plot(time_cc,sT_plusConverg_minusPairwise_cc_mean(a,:),'k--')
%     plot(time_cc,sT_minusConverg_minusPairwise_cc_mean(a,:),'g--')
%     legend('+c+p','-c+p','+c-p','-c-p')
%     xlim([-1 1])
%     ylim([-.5 1])
% end
% 
% pause
% close

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
%     hold on
%     plot([1:4],sT_corrCoefs(a,:),'-+')
%     
%     legend('sn','st')
%     set(gca,'XTick',[1:4],'XTickLabel',{'+c+p','-c+p','+c-p','-c-p'})
%     xlabel('input correlations')
%     ylabel('output sn correlation coef')
% end
% 
% pause
% close

figure % bubble plots
for a = 1:numBars ;
    plot(sn_plusConverg_plusPairwise_ellipseX(a,:),sn_plusConverg_plusPairwise_ellipseY(a,:),'k')
    hold on
    plot(sn_minusConverg_plusPairwise_ellipseX(a,:),sn_minusConverg_plusPairwise_ellipseY(a,:),'g')
    plot(sn_plusConverg_minusPairwise_ellipseX(a,:),sn_plusConverg_minusPairwise_ellipseY(a,:),'c')
    plot(sn_minusConverg_minusPairwise_ellipseX(a,:),sn_minusConverg_minusPairwise_ellipseY(a,:),'y')
end
legend('+c+p','-c+p','+c-p','-c-p')

% FOR IGOR

% % spike number and variance
% identifier = ['sn1SimMean','id','cell',num2str(A)] ;
% ForIgor.(identifier) = sn1_sim_mean ; 
% 
% identifier = ['sn2SimMean','id','cell',num2str(A)] ;
% ForIgor.(identifier) = sn2_sim_mean ; 
% 
% identifier = ['sn1SimVar','id','cell',num2str(A)] ;
% ForIgor.(identifier) = sn1_sim_var ; 
% 
% identifier = ['sn2SimVar','id','cell',num2str(A)] ;
% ForIgor.(identifier) = sn2_sim_var;

% spike number and variance per bar
% for a=1:numBars
%     identifier = ['sn1SimMean',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = sn1_sim_mean(a) ; 
% 
%     identifier = ['sn2SimMean',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = sn2_sim_mean(a) ; 
% 
%     identifier = ['sn1SimVar',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = sn1_sim_var(a) ; 
% 
%     identifier = ['sn2SimVar',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = sn2_sim_var(a);
% 
%     identifier = ['sn1ShuffMean',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = sn1_shuff_mean(a) ; 
% 
%     identifier = ['sn2ShuffMean',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = sn2_shuff_mean(a) ; 
% 
%     identifier = ['sn1ShuffVar',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = sn1_shuff_var(a) ; 
% 
%     identifier = ['sn2ShuffVar',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = sn2_shuff_var(a);
% 
% end
% 
% % spike correlations
% for a = 1:numBars ;
%     identifier = ['snCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = sn_corrCoefs(a,:);
%     
%     identifier = ['stCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = sT_corrCoefs(a,:);
%     
%     identifier = ['snCorrsStd',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = sn_corrCoefsStd(a,:);
% end

% make room in memmory
clearvars -except ForIgor Parameters id A excg* inhg* sT_corrCoefs sn_corrCoefs num* prePnts barPnts minusPairwiseShift SI time Voltage_mean SpikeThreshold_mean data SpikePnts absRef InhRev

%% spike generation and correlations through an LIF or EIF model
%
samplerate = 1/SI ;
params.Eexc = 0 ;          % reversal potential for excitatory current
params.Einh = -60;
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

numberTrials_optimized = 60 ;
for b=1:numberTrials_optimized ;
    % optimize Iadd to capture the voltage mean accurately 
    for a=1:length(Iadd_range) ;
        LIFv = LIFmodelGplusI(excg(b,:),inhg(b,:),0,Iadd_range(a),samplerate,params) ; % generate voltage data with LIF model
        MeanVerror(a) = (mean(LIFv) - mean(data(b,:)))^2 ; 
    end
    [tempV,tempi] = min(MeanVerror) ;
    Iadd_opt(b) = Iadd_range(tempi) ;

    % optimize Gleak to capture the voltage variance accurately
    for a=1:length(Gleak_range) ;
        LIFv = LIFmodelGplusI(excg(b,:),inhg(b,:),Gleak_range(a),Iadd_opt(b),samplerate,params) ; % generate voltage data with LIF model
        varVerror(a) = (var(LIFv) - var(data(b,:)))^2 ; 
    end
    [tempV,tempi] = min(varVerror) ;
    Gleak_opt(b) = Gleak_range(tempi) ;

    % optimize spike Thresh to capture the spike number accuratley
    for a =1:length(Vthresh_range) ;
        params.Vthresh = Vthresh_range(a) ;
        LIFv = LIFmodelGplusI(excg(b,:),inhg(b,:),Gleak_opt(b),Iadd_opt(b),samplerate,params) ; 
        snerror(a) = (length(find(LIFv==50)) - length(SpikePnts{b}))^2 ; 
    end
    [tempV,tempi] = min(snerror) ;
    Vthresh_opt(b) = Vthresh_range(tempi) ;
end
params.Vthresh = mean(Vthresh_opt) ;

% run LIF model with optimized parameters
LIFv = LIFmodelGplusI(excg,inhg,mean(Gleak_opt),mean(Iadd_opt),samplerate,params) ;

% % EIF model
% SafeRef = .2 ; % seconds post spike that voltage is still impacted by post spike intrinsic g
% SafeRefPnts = round(SafeRef/SI(1)) ; 
% 
% totalNumSpikes = 0 ;
% for a=1:length(SpikePnts) ; % for every spike epoch
%     totalNumSpikes = totalNumSpikes + length(SpikePnts{a}) ;
% end
% preSpikeData = nans(totalNumSpikes,SafeRefPnts) ; % preallocate for speed
% preSpikeGexc = nans(totalNumSpikes,SafeRefPnts) ; % preallocate for speed
% preSpikeGinh = nans(totalNumSpikes,SafeRefPnts) ; % preallocate for speed
% 
% r=0 ;
% for a=1:length(SpikePnts) ; % for every spike epoch
%     if ~isempty(SpikePnts{a}) ; % if there are spikes
%         for b=1:length(SpikePnts{a}) ; % for every spike
%             if (b==1 && SpikePnts{a}(b)>SafeRefPnts) || (b>1 && (SpikePnts{a}(b)-SpikePnts{a}(b-1))>SafeRefPnts) % if your spike is not preceded by another spike within the safe refractory
%                 r=r+1 ;
%                 preSpikeData(r,:) = data(a,SpikePnts{a}(b)+1-SafeRefPnts:SpikePnts{a}(b)) ;
%                 preSpikeGexc(r,:) = excg(a,SpikePnts{a}(b)+1-SafeRefPnts:SpikePnts{a}(b)) ;
%                 preSpikeGinh(r,:) = inhg(a,SpikePnts{a}(b)+1-SafeRefPnts:SpikePnts{a}(b)) ;
%             end
%         end
%     end
% end
% preSpikeData = preSpikeData(1:r,:) ;
% preSpikeGexc = preSpikeGexc(1:r,:)  ;
% preSpikeGinh = preSpikeGinh(1:r,:)  ; 
% 
% [Voltage_unique,fV,fV_sem] = fVgClampFinder(preSpikeData,preSpikeGexc,preSpikeGinh,params.Eexc,params.Einh,0, params.cap,2,SI(1)) ;
% 
% LIFv = EIFmodelGplusGleakAndIadd(excg,inhg,0,-100,Voltage_unique,fV,samplerate,params) ;
% LIFv = EIFmodelG(excg,inhg,Voltage_unique,fV,samplerate,params) ;

% mean and variance of spike number in data and model
for a=1:size(data,1) ;
    snData(a) = length(SpikePnts{a}) ;
    snLIF(a) = length(find(LIFv(a,:)==50)) ;
end
snData_mean = mean(snData) ;
snLIF_mean = mean(snLIF) ;
snData_var = var(snData) ;
snLIF_var = var(snLIF) ;


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
numBars = 3 ; % number of bars presented while recording conductances
prePnts = 2000 ; % pnts to avoid as g starts to be injected
barPnts = floor((length(data)-prePnts)/numBars) ; 

for a=1:numBars ; % preallocate variables to improve memory usage and speed
    spikeTrain1_sim_bar{a} = nans(numTrials,barPnts) ;
    spikeTrain2_sim_bar{a} = nans(numTrials,barPnts) ;
    spikeTrain1_shuff_bar{a} = nans(numTrials,barPnts) ;
    spikeTrain2_shuff_bar{a} = nans(numTrials,barPnts) ;
end

for a=1:numBars ;
    spikeTrain1_sim_bar{a} = spikeTrain1_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    spikeTrain2_sim_bar{a} = spikeTrain2_sim(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    spikeTrain1_shuff_bar{a} = spikeTrain1_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    spikeTrain2_shuff_bar{a} = spikeTrain2_shuff(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
end

% spike number for each bar
sn1_sim = nans(numTrials,numBars) ; % preallocate variables to improve memory usage and speed
sn2_sim = nans(numTrials,numBars) ;
sn1_shuff = nans(numTrials,numBars) ;
sn2_shuff = nans(numTrials,numBars) ;
sn1_sim_mean = nans(1,numBars) ;
sn2_sim_mean = nans(1,numBars) ;
sn1_shuff_mean = nans(1,numBars) ;
sn2_shuff_mean = nans(1,numBars) ;

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

% smoothed spike trains
for a=1:numBars ; % preallocate variables to improve memory usage and speed
    spikeTrain1_smth_sim_bar{a} = nans(numTrials,barPnts) ;
    spikeTrain2_smth_sim_bar{a} = nans(numTrials,barPnts) ;
    spikeTrain1_smth_shuff_bar{a} = nans(numTrials,barPnts) ;
    spikeTrain2_smth_shuff_bar{a} = nans(numTrials,barPnts) ;    
end

smthFctr = .050 ; % sec
smthFctr = round(smthFctr/SI) ;
for a=1:numBars ;
    for b=1:numTrials ;
        spikeTrain1_smth_sim_bar{a}(b,:) = smooth(spikeTrain1_sim_bar{a}(b,:),smthFctr) ;
        spikeTrain2_smth_sim_bar{a}(b,:) = smooth(spikeTrain2_sim_bar{a}(b,:),smthFctr) ;
        spikeTrain1_smth_shuff_bar{a}(b,:) = smooth(spikeTrain1_shuff_bar{a}(b,:),smthFctr) ;
        spikeTrain2_smth_shuff_bar{a}(b,:) = smooth(spikeTrain2_shuff_bar{a}(b,:),smthFctr) ;    
    end
end

% spike number residuals (local correction)
sn1_sim_res = nans(numTrials-2,numBars) ; % preallocate variables to improve memory usage and speed
sn2_sim_res = nans(numTrials-2,numBars) ;
sn1_shuff_res = nans(numTrials-2,numBars) ;
sn2_shuff_res = nans(numTrials-2,numBars) ;

for a = 2:numTrials-1 ;
    sn1_sim_res(a-1,:) = sn1_sim(a,:) - (sn1_sim(a-1,:)+sn1_sim(a+1,:))/2 ;
    sn2_sim_res(a-1,:) = sn2_sim(a,:) - (sn2_sim(a-1,:)+sn2_sim(a+1,:))/2 ;
    sn1_shuff_res(a-1,:) = sn1_shuff(a,:) - (sn1_shuff(a-1,:)+sn1_shuff(a+1,:))/2 ;
    sn2_shuff_res(a-1,:) = sn2_shuff(a,:) - (sn2_shuff(a-1,:)+sn2_shuff(a+1,:))/2 ;
end

% smoothed spike train residuals
spikeTrain1_smth_sim_bar_res = cell(1,numBars) ; % preallocate variables to improve memory usage and speed
spikeTrain2_smth_sim_bar_res = cell(1,numBars) ;
spikeTrain1_smth_shuff_bar_res = cell(1,numBars) ;
spikeTrain2_smth_shuff_bar_res = cell(1,numBars) ;    

for a=1:numBars ;
    for b=2:numTrials-1  ;
        spikeTrain1_smth_sim_bar_res{a}(b-1,:) = spikeTrain1_smth_sim_bar{a}(b,:) - (spikeTrain1_smth_sim_bar{a}(b-1,:) + spikeTrain1_smth_sim_bar{a}(b+1,:))./2 ;
        spikeTrain2_smth_sim_bar_res{a}(b-1,:) = spikeTrain2_smth_sim_bar{a}(b,:) - (spikeTrain2_smth_sim_bar{a}(b-1,:) + spikeTrain2_smth_sim_bar{a}(b+1,:))./2 ;
        spikeTrain1_smth_shuff_bar_res{a}(b-1,:) = spikeTrain1_smth_shuff_bar{a}(b,:) - (spikeTrain1_smth_shuff_bar{a}(b-1,:) + spikeTrain1_smth_shuff_bar{a}(b+1,:))./2 ;
        spikeTrain2_smth_shuff_bar_res{a}(b-1,:) = spikeTrain2_smth_shuff_bar{a}(b,:) - (spikeTrain2_smth_shuff_bar{a}(b-1,:) + spikeTrain2_smth_shuff_bar{a}(b+1,:))./2 ;    
    end
end

% sn unbiased variance of spike number from residuals
sn1_sim_var = sum(sn1_sim_res.^2,1)./(size(sn1_sim_res,1)-1) ;
sn2_sim_var = sum(sn2_sim_res.^2,1)./(size(sn2_sim_res,1)-1) ;
sn1_shuff_var = sum(sn1_shuff_res.^2,1)./(size(sn1_shuff_res,1)-1) ;
sn2_shuff_var = sum(sn2_shuff_res.^2,1)./(size(sn2_shuff_res,1)-1) ;

% spike number correlation (linear) for each bar
lincoefsigMin = .05 ;
minusPairwiseShift = 3 ; % number of trials to shift to get rid of pairwise correlations
for a=1:numBars ;
    [Corr,p] =  corrcoef(sn1_sim_res(:,a),sn2_sim_res(:,a)) ; % + converging + pairwise correlations 
    sn_plusConverg_plusPairwise_lincoef(a) = Corr(1,2) ;
    sn_plusConverg_plusPairwise_lincoefsig(a) = p(1,2) ;
%     if sn_plusConverg_plusPairwise_lincoefsig(a)>lincoefsigMin ;
%         sn_plusConverg_plusPairwise_lincoef(a) = 0 ;
%     end

    [Corr,p] =  corrcoef(sn1_shuff_res(:,a),sn2_shuff_res(:,a)) ; % - converging + pairwise correlations
    sn_minusConverg_plusPairwise_lincoef(a) = Corr(1,2) ;
    sn_minusConverg_plusPairwise_lincoefsig(a) = p(1,2) ;
%     if sn_minusConverg_plusPairwise_lincoefsig(a)>lincoefsigMin ;
%         sn_minusConverg_plusPairwise_lincoef(a) = 0 ;
%     end
    
    [Corr,p] =  corrcoef(sn1_sim_res(:,a),circshift(sn2_sim_res(:,a),[minusPairwiseShift,0])) ; % + converging - pairwise correlations 
    sn_plusConverg_minusPairwise_lincoef(a) = Corr(1,2) ;
    sn_plusConverg_minusPairwise_lincoefsig(a) = p(1,2) ;
%     if sn_plusConverg_minusPairwise_lincoefsig(a)>lincoefsigMin ;
%         sn_plusConverg_minusPairwise_lincoef(a) = 0 ;
%     end
    
    [Corr,p] =  corrcoef(sn1_shuff_res(:,a),circshift(sn2_shuff_res(:,a),[minusPairwiseShift,0])) ; % - converging - pairwise correlations 
    sn_minusConverg_minusPairwise_lincoef(a) = Corr(1,2) ;
    sn_minusConverg_minusPairwise_lincoefsig(a) = p(1,2) ;
%     if sn_minusConverg_minusPairwise_lincoefsig(a)>lincoefsigMin ;
%         sn_minusConverg_minusPairwise_lincoef(a) = 0 ;
%     end

    LIFsn_corrCoefs(a,:) = [sn_plusConverg_plusPairwise_lincoef(a),sn_minusConverg_plusPairwise_lincoef(a),...
    sn_plusConverg_minusPairwise_lincoef(a),sn_minusConverg_minusPairwise_lincoef(a)] ;
end

% smoothed spike train correlations
for a=1:numBars ; % preallocate variables to improve memory usage and speed
    spikeTrain2_smth_sim_bar_res_circshift{a} = nans(numTrials-2,barPnts) ; 
    spikeTrain2_smth_shuff_bar_res_circshift{a} = nans(numTrials-2,barPnts) ;        
    sT_plusConverg_plusPairwise_cc{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_minusConverg_plusPairwise_cc{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_plusConverg_minusPairwise_cc{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_minusConverg_minusPairwise_cc{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_plusConverg_plusPairwise_c{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_minusConverg_plusPairwise_c{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_plusConverg_minusPairwise_c{a} = nans(numTrials-2,barPnts*2-1) ;
    sT_minusConverg_minusPairwise_c{a} = nans(numTrials-2,barPnts*2-1) ;
end
sT_plusConverg_plusPairwise_ccPeak = nans(numBars,numTrials-2) ;
sT_minusConverg_plusPairwise_ccPeak = nans(numBars,numTrials-2) ;
sT_plusConverg_minusPairwise_ccPeak = nans(numBars,numTrials-2) ;
sT_minusConverg_minusPairwise_ccPeak = nans(numBars,numTrials-2) ;
sT_plusConverg_plusPairwise_cPeak = nans(numBars,numTrials-2) ;
sT_minusConverg_plusPairwise_cPeak = nans(numBars,numTrials-2) ;
sT_plusConverg_minusPairwise_cPeak = nans(numBars,numTrials-2) ;
sT_minusConverg_minusPairwise_cPeak = nans(numBars,numTrials-2) ;
sT_plusConverg_plusPairwise_cc_mean = nans(numTrials-2,barPnts*2-1) ; 
sT_minusConverg_plusPairwise_cc_mean = nans(numTrials-2,barPnts*2-1) ;
sT_plusConverg_minusPairwise_cc_mean = nans(numTrials-2,barPnts*2-1) ;
sT_minusConverg_minusPairwise_cc_mean = nans(numTrials-2,barPnts*2-1) ;
sT_plusConverg_plusPairwise_c_mean = nans(numTrials-2,barPnts*2-1) ; 
sT_minusConverg_plusPairwise_c_mean = nans(numTrials-2,barPnts*2-1) ;
sT_plusConverg_minusPairwise_c_mean = nans(numTrials-2,barPnts*2-1) ;
sT_minusConverg_minusPairwise_c_mean = nans(numTrials-2,barPnts*2-1) ;
sT_plusConverg_plusPairwise_cc_mean_Peak = nans(1,numBars) ; 
sT_minusConverg_plusPairwise_cc_mean_Peak = nans(1,numBars) ;
sT_plusConverg_minusPairwise_cc_mean_Peak = nans(1,numBars) ;
sT_minusConverg_minusPairwise_cc_mean_Peak = nans(1,numBars) ;    
sT_plusConverg_plusPairwise_c_mean_Peak = nans(1,numBars) ; 
sT_minusConverg_plusPairwise_c_mean_Peak = nans(1,numBars) ;
sT_plusConverg_minusPairwise_c_mean_Peak = nans(1,numBars) ;
sT_minusConverg_minusPairwise_c_mean_Peak = nans(1,numBars) ; 


for a=1:numBars ;
    spikeTrain2_smth_sim_bar_res_circshift{a} = circshift(spikeTrain2_smth_sim_bar_res{a},[minusPairwiseShift,0]) ;
    spikeTrain2_smth_shuff_bar_res_circshift{a} = circshift(spikeTrain2_smth_shuff_bar_res{a},[minusPairwiseShift,0]) ;
    
    for b=1:numTrials-2  ;
        % corr coef
        sT_plusConverg_plusPairwise_cc{a}(b,:) = xcorr(spikeTrain1_smth_sim_bar_res{a}(b,:),spikeTrain2_smth_sim_bar_res{a}(b,:),'coeff') ;
        sT_minusConverg_plusPairwise_cc{a}(b,:) = xcorr(spikeTrain1_smth_shuff_bar_res{a}(b,:),spikeTrain2_smth_shuff_bar_res{a}(b,:),'coeff') ;
        sT_plusConverg_minusPairwise_cc{a}(b,:) = xcorr(spikeTrain1_smth_sim_bar_res{a}(b,:),spikeTrain2_smth_sim_bar_res_circshift{a}(b,:),'coeff') ;
        sT_minusConverg_minusPairwise_cc{a}(b,:) = xcorr(spikeTrain1_smth_shuff_bar_res{a}(b,:),spikeTrain2_smth_shuff_bar_res_circshift{a}(b,:),'coeff') ;
    
        % corr normalized by length only
        sT_plusConverg_plusPairwise_c{a}(b,:) = xcorr(spikeTrain1_smth_sim_bar_res{a}(b,:),spikeTrain2_smth_sim_bar_res{a}(b,:),'unbiased') ;
        sT_minusConverg_plusPairwise_c{a}(b,:) = xcorr(spikeTrain1_smth_shuff_bar_res{a}(b,:),spikeTrain2_smth_shuff_bar_res{a}(b,:),'unbiased') ;
        sT_plusConverg_minusPairwise_c{a}(b,:) = xcorr(spikeTrain1_smth_sim_bar_res{a}(b,:),spikeTrain2_smth_sim_bar_res_circshift{a}(b,:),'unbiased') ;
        sT_minusConverg_minusPairwise_c{a}(b,:) = xcorr(spikeTrain1_smth_shuff_bar_res{a}(b,:),spikeTrain2_smth_shuff_bar_res_circshift{a}(b,:),'unbiased') ;
        
        % peaks of cc
        sT_plusConverg_plusPairwise_ccPeak(a,b) = CCpeakFinder(sT_plusConverg_plusPairwise_cc{a}(b,:)) ;
        sT_minusConverg_plusPairwise_ccPeak(a,b) = CCpeakFinder(sT_minusConverg_plusPairwise_cc{a}(b,:)) ;
        sT_plusConverg_minusPairwise_ccPeak(a,b) = CCpeakFinder(sT_plusConverg_minusPairwise_cc{a}(b,:)) ;
        sT_minusConverg_minusPairwise_ccPeak(a,b) = CCpeakFinder(sT_minusConverg_minusPairwise_cc{a}(b,:)) ;
        
        % peaks of c
        sT_plusConverg_plusPairwise_cPeak(a,b) = CCpeakFinder(sT_plusConverg_plusPairwise_c{a}(b,:)) ;
        sT_minusConverg_plusPairwise_cPeak(a,b) = CCpeakFinder(sT_minusConverg_plusPairwise_c{a}(b,:)) ;
        sT_plusConverg_minusPairwise_cPeak(a,b) = CCpeakFinder(sT_plusConverg_minusPairwise_c{a}(b,:)) ;
        sT_minusConverg_minusPairwise_cPeak(a,b) = CCpeakFinder(sT_minusConverg_minusPairwise_c{a}(b,:)) ;
    end
    % mean corr coef
    sT_plusConverg_plusPairwise_cc_mean(a,:) = nanmean(sT_plusConverg_plusPairwise_cc{a}) ; % trials where residual = 0 have no cc (no residual = no noise = no ability to correlate noise)
    sT_minusConverg_plusPairwise_cc_mean(a,:) = nanmean(sT_minusConverg_plusPairwise_cc{a}) ;
    sT_plusConverg_minusPairwise_cc_mean(a,:) = nanmean(sT_plusConverg_minusPairwise_cc{a}) ;
    sT_minusConverg_minusPairwise_cc_mean(a,:) = nanmean(sT_minusConverg_minusPairwise_cc{a}) ;
    
    % mean corr 
    sT_plusConverg_plusPairwise_c_mean(a,:) = nanmean(sT_plusConverg_plusPairwise_c{a}) ; 
    sT_minusConverg_plusPairwise_c_mean(a,:) = nanmean(sT_minusConverg_plusPairwise_c{a}) ;
    sT_plusConverg_minusPairwise_c_mean(a,:) = nanmean(sT_plusConverg_minusPairwise_c{a}) ;
    sT_minusConverg_minusPairwise_c_mean(a,:) = nanmean(sT_minusConverg_minusPairwise_c{a}) ;
    
    % peak of mean corr coef
    sT_plusConverg_plusPairwise_cc_mean_Peak(a) = CCpeakFinder(sT_plusConverg_plusPairwise_cc_mean(a,:)) ; 
    sT_minusConverg_plusPairwise_cc_mean_Peak(a) = CCpeakFinder(sT_minusConverg_plusPairwise_cc_mean(a,:)) ;
    sT_plusConverg_minusPairwise_cc_mean_Peak(a) = CCpeakFinder(sT_plusConverg_minusPairwise_cc_mean(a,:)) ;
    sT_minusConverg_minusPairwise_cc_mean_Peak(a) = CCpeakFinder(sT_minusConverg_minusPairwise_cc_mean(a,:)) ;    
    
    % peak of mean corr
    sT_plusConverg_plusPairwise_c_mean_Peak(a) = CCpeakFinder(sT_plusConverg_plusPairwise_c_mean(a,:)) ; 
    sT_minusConverg_plusPairwise_c_mean_Peak(a) = CCpeakFinder(sT_minusConverg_plusPairwise_c_mean(a,:)) ;
    sT_plusConverg_minusPairwise_c_mean_Peak(a) = CCpeakFinder(sT_plusConverg_minusPairwise_c_mean(a,:)) ;
    sT_minusConverg_minusPairwise_c_mean_Peak(a) = CCpeakFinder(sT_minusConverg_minusPairwise_c_mean(a,:)) ;   
    
    % concatinated corr coefs
    LIFsT_corrCoefs(a,:) = [sT_plusConverg_plusPairwise_cc_mean_Peak(a),sT_minusConverg_plusPairwise_cc_mean_Peak(a),...
    sT_plusConverg_minusPairwise_cc_mean_Peak(a),sT_minusConverg_minusPairwise_cc_mean_Peak(a)] ;

    % concatinated corrs
    LIFsT_corr(a,:) = [sT_plusConverg_plusPairwise_c_mean_Peak(a),sT_minusConverg_plusPairwise_c_mean_Peak(a),...
    sT_plusConverg_minusPairwise_c_mean_Peak(a),sT_minusConverg_minusPairwise_c_mean_Peak(a)] ;
end


% figure % comparing dynamic clamp and model voltage predictions
% for a=1:size(data,1) ;
%     plot(LIFv(a,:),'r')
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
% 
%     plot(SpikePnts{a},a,'b*')
%     hold on
%     plot(find(LIFv(a,:)==50),a,'ro')
% end

figure % comparing dynamic clamp and model spike number mean and variance
plot(snData,'b')
hold on
plot(snLIF,'r')

% figure % quantitative impact of converging and pairwise correlations on spike correlations
% for a = 1:numBars ;
%     subplot(1,numBars,a)
%     plot([1:4],LIFsn_corrCoefs(a,:),'-o')
%     hold on
%     plot([1:4],LIFsT_corrCoefs(a,:),'-+')
%     
%     legend('sn','st')
%     set(gca,'XTick',[1:4],'XTickLabel',{'+c+p','-c+p','+c-p','-c-p'})
%     xlabel('input correlations')
%     ylabel('output sn correlation coef')
% end
% 
% pause
% close

% FOR IGOR

% LIF for all conditions
% for a = 1:numBars ;
%     identifier = ['LIFsnCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = LIFsn_corrCoefs(a,:);
%     
%     identifier = ['LIFstCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = LIFsT_corrCoefs(a,:);    
% end

% LIF for +c+p vs -c+p only
identifier = ['snCorrsActual',id,'cell',num2str(A)] ;
ForIgor.(identifier) = sn_corrCoefs(1:numBars*2) ;

identifier = ['snCorrsLIFPred',id,'cell',num2str(A)] ;
ForIgor.(identifier) = LIFsn_corrCoefs(1:numBars*2);
    

%
% make room in memmory
clearvars -except ForIgor Parameters id A excg* inhg* sT_corrCoefs sn_corrCoefs num* prePnts barPnts minusPairwiseShift SI time Voltage_mean InhRev

%% noise correlations and variance in presented conductances  
%{

% conductance residuals (local correction)
for b=1:numBars ; % preallocate variables to improve memory usage and speed 
    excg1_sim_bar_res{b} = nans(numTrials-2,barPnts) ;
    excg2_sim_bar_res{b} = nans(numTrials-2,barPnts) ;
    excg1_shuff_bar_res{b} = nans(numTrials-2,barPnts) ;
    excg2_shuff_bar_res{b} = nans(numTrials-2,barPnts) ;
    inhg1_sim_bar_res{b} = nans(numTrials-2,barPnts) ;
    inhg2_sim_bar_res{b} = nans(numTrials-2,barPnts) ;
    inhg1_shuff_bar_res{b} = nans(numTrials-2,barPnts) ;
    inhg2_shuff_bar_res{b} = nans(numTrials-2,barPnts) ;
end

for b=1:numBars ;
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
    c1_plusConverge{b} = nans(numTrials-2,barPnts*2-1) ;
    c2_plusConverge{b} = nans(numTrials-2,barPnts*2-1) ;
    c1_minusConverge{b} = nans(numTrials-2,barPnts*2-1) ;
    c2_minusConverge{b} = nans(numTrials-2,barPnts*2-1) ;
    cExc_plusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    cInh_plusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;  
    cExc_minusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    cInh_minusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;   
    cExc_plusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    cInh_plusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    cExc_minusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    cInh_minusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ; 
    cExcInh_plusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    cInhExc_plusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    cExcInh_minusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    cInhExc_minusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    cExcInh_plusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    cInhExc_plusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    cExcInh_minusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    cInhExc_minusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ; 
    excg1_sim_bar_res_ac{b} = nans(numTrials-2,barPnts*2-1) ;
    inhg1_sim_bar_res_ac{b} = nans(numTrials-2,barPnts*2-1) ;
    excg2_sim_bar_res_ac{b} = nans(numTrials-2,barPnts*2-1) ;
    inhg2_sim_bar_res_ac{b} = nans(numTrials-2,barPnts*2-1) ;
    excg1_shuff_bar_res_ac{b} = nans(numTrials-2,barPnts*2-1) ;
    inhg1_shuff_bar_res_ac{b} = nans(numTrials-2,barPnts*2-1) ;
    excg2_shuff_bar_res_ac{b} = nans(numTrials-2,barPnts*2-1) ;
    inhg2_shuff_bar_res_ac{b} = nans(numTrials-2,barPnts*2-1) ;
    excg2_sim_bar_res_circshift_ac{b} = nans(numTrials-2,barPnts*2-1) ;
    inhg2_sim_bar_res_circshift_ac{b} = nans(numTrials-2,barPnts*2-1) ;
    excg2_shuff_bar_res_circshift_ac{b} = nans(numTrials-2,barPnts*2-1) ;
    inhg2_shuff_bar_res_circshift_ac{b} = nans(numTrials-2,barPnts*2-1) ;
end
c1_plusConverge_peak = nans(numBars,numTrials-2) ; 
c2_plusConverge_peak = nans(numBars,numTrials-2) ;  
c1_minusConverge_peak = nans(numBars,numTrials-2) ; 
c2_minusConverge_peak = nans(numBars,numTrials-2) ; 
cExc_plusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ; 
cInh_plusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ;         
cExc_minusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ; 
cInh_minusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ;        
cExc_plusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ; 
cInh_plusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ;        
cExc_minusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ; 
cInh_minusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ;         
cExcInh_plusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ; 
cInhExc_plusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ;  
cExcInh_minusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ; 
cInhExc_minusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ;       
cExcInh_plusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ; 
cInhExc_plusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ; 
cExcInh_minusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ; 
cInhExc_minusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ;  
cSpikePred_plusConverg_plusPairwise = nans(numBars,numTrials-2) ; 
cSpikePred_minusConverg_plusPairwise = nans(numBars,numTrials-2) ; 
cSpikePred_plusConverg_minusPairwise = nans(numBars,numTrials-2) ; 
cSpikePred_minusConverg_minusPairwise = nans(numBars,numTrials-2) ;         
c1_plusConverge_mean = nans(numBars,barPnts*2-1) ; 
c2_plusConverge_mean = nans(numBars,barPnts*2-1) ;    
c1_minusConverge_mean = nans(numBars,barPnts*2-1) ; 
c2_minusConverge_mean = nans(numBars,barPnts*2-1) ;     
cExc_plusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cInh_plusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cExc_minusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cInh_minusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cExc_plusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cInh_plusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cExc_minusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cInh_minusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ;     
cExcInh_plusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cInhExc_plusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cExcInh_minusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cInhExc_minusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cExcInh_plusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cInhExc_plusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cExcInh_minusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ; 
cInhExc_minusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ;  


for b=1:numBars ;
    excg2_sim_bar_res_circshift{b} = circshift(excg2_sim_bar_res{b},[minusPairwiseShift,0]) ;
    inhg2_sim_bar_res_circshift{b} = circshift(inhg2_sim_bar_res{b},[minusPairwiseShift,0]) ;
    excg2_shuff_bar_res_circshift{b} = circshift(excg2_shuff_bar_res{b},[minusPairwiseShift,0]) ;
    inhg2_shuff_bar_res_circshift{b} = circshift(inhg2_shuff_bar_res{b},[minusPairwiseShift,0]) ;
    
    for a = 1:numTrials-2 ; % for each residual     
        
        % converging correlations unormalized
        c1_plusConverge{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),inhg1_sim_bar_res{b}(a,:),'unbiased') ;
        c2_plusConverge{b}(a,:) = xcorr(excg2_sim_bar_res{b}(a,:),inhg2_sim_bar_res{b}(a,:),'unbiased') ;
        
        c1_minusConverge{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),inhg1_shuff_bar_res{b}(a,:),'unbiased') ;
        c2_minusConverge{b}(a,:) = xcorr(excg2_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res{b}(a,:),'unbiased') ;
    

        % pairwise ee,ii correlations unormalized
        cExc_plusConverg_plusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),excg2_sim_bar_res{b}(a,:),'unbiased') ; % + converging + pairwise
        cInh_plusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res{b}(a,:),'unbiased') ;
        
        cExc_minusConverg_plusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res{b}(a,:),'unbiased') ; % -converging + pairwise
        cInh_minusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res{b}(a,:),'unbiased') ;
        
        cExc_plusConverg_minusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),excg2_sim_bar_res_circshift{b}(a,:),'unbiased') ; % + converging - pairwise
        cInh_plusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res_circshift{b}(a,:),'unbiased') ;
        
        cExc_minusConverg_minusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res_circshift{b}(a,:),'unbiased') ; % -converging - pairwise
        cInh_minusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res_circshift{b}(a,:),'unbiased') ; 
  
        % pairwise ei,ie correlations unormalized
        cExcInh_plusConverg_plusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res{b}(a,:),'unbiased') ; % + converging + pairwise
        cInhExc_plusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),excg2_sim_bar_res{b}(a,:),'unbiased') ;
        
        cExcInh_minusConverg_plusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res{b}(a,:),'unbiased') ; % -converging + pairwise
        cInhExc_minusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res{b}(a,:),'unbiased') ;
        
        cExcInh_plusConverg_minusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res_circshift{b}(a,:),'unbiased') ; % + converging - pairwise
        cInhExc_plusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),excg2_sim_bar_res_circshift{b}(a,:),'unbiased') ;
        
        cExcInh_minusConverg_minusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res_circshift{b}(a,:),'unbiased') ; % -converging - pairwise
        cInhExc_minusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res_circshift{b}(a,:),'unbiased') ;    
        
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
        excg1_sim_bar_res_ac{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),'unbiased') ;
        inhg1_sim_bar_res_ac{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),'unbiased') ;
        excg2_sim_bar_res_ac{b}(a,:) = xcorr(excg2_sim_bar_res{b}(a,:),'unbiased') ;
        inhg2_sim_bar_res_ac{b}(a,:) = xcorr(inhg2_sim_bar_res{b}(a,:),'unbiased') ;
        excg1_shuff_bar_res_ac{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),'unbiased') ;
        inhg1_shuff_bar_res_ac{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),'unbiased') ;
        excg2_shuff_bar_res_ac{b}(a,:) = xcorr(excg2_shuff_bar_res{b}(a,:),'unbiased') ;
        inhg2_shuff_bar_res_ac{b}(a,:) = xcorr(inhg2_shuff_bar_res{b}(a,:),'unbiased') ;
        excg2_sim_bar_res_circshift_ac{b}(a,:) = xcorr(excg2_sim_bar_res_circshift{b}(a,:),'unbiased') ;
        inhg2_sim_bar_res_circshift_ac{b}(a,:) = xcorr(inhg2_sim_bar_res_circshift{b}(a,:),'unbiased') ;
        excg2_shuff_bar_res_circshift_ac{b}(a,:) = xcorr(excg2_shuff_bar_res_circshift{b}(a,:),'unbiased') ;
        inhg2_shuff_bar_res_circshift_ac{b}(a,:) = xcorr(inhg2_shuff_bar_res_circshift{b}(a,:),'unbiased') ;
        
 
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

time_cc = SI*([1:length(c1_plusConverge_mean)]-(length(c1_plusConverge_mean)+1)/2) ;

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

% solving for free parameter
%{
for a=1:numBars; % for every bar
    for b=1:4 ; % for every set pairings   
    [alpha_sn(a,b),fval_sn(a,b),ef_sn(a,b)] = fsolve(@(alpha) PairCorrAnalyticalSolutionfit(alpha,Ve1e2(a,b),Vi1i2(a,b),Ve1i2(a,b),Vi1e2(a,b),Ve1(a,b),Vi1(a,b),Ve1i1(a,b),Ve2(a,b),Vi2(a,b),Ve2i2(a,b),sn_corrCoefs(a,b)),0.5) ;
    [alpha_st(a,b),fval_st(a,b),ef_st(a,b)] = fsolve(@(alpha) PairCorrAnalyticalSolutionfit(alpha,Ve1e2(a,b),Vi1i2(a,b),Ve1i2(a,b),Vi1e2(a,b),Ve1(a,b),Vi1(a,b),Ve1i1(a,b),Ve2(a,b),Vi2(a,b),Ve2i2(a,b),sT_corrCoefs(a,b)),0.5) ;
    end
end

% fraction of g sets that were solvable with linear equation
FracFoundAlpha_sn = mean(ef_sn(:)==1) ;
FracFoundAlpha_st = mean(ef_st(:)==1) ;

% alpha values found when solvable
alphaHistX = [-1:.2:2] ;

i1 = find(ef_sn==1) ;
alpha_snHist = hist(alpha_sn(i1),alphaHistX) ;

i1 = find(ef_st==1) ;
alpha_stHist = hist(alpha_st(i1),alphaHistX) ;

%}

% assume an alpha and predict Pss
numConditions = 4 ; 
alpha_range = [0:.05:1] ;
for iter = 1:length(alpha_range) ;
    alpha = alpha_range(iter) ;
    for a=1:numBars; % for every bar
        for b=1:numConditions ; % for every set pairings   
            numerator = Ve1e2(a,b) + alpha^2*Vi1i2(a,b) - alpha*Ve1i2(a,b) - alpha*Vi1e2(a,b) ;
            denominator = sqrt((Ve1(a,b) + alpha^2*Vi1(a,b) - 2*alpha*Ve1i1(a,b))*(Ve2(a,b) + alpha^2*Vi2(a,b) - 2*alpha*Ve2i2(a,b))) ;
            Pss_predicted(a,b) = numerator/denominator ;
            Pss_predictedNoCross(a,b) = (Ve1e2(a,b) + alpha^2*Vi1i2(a,b))/sqrt((Ve1(a,b) + alpha^2*Vi1(a,b))*(Ve2(a,b) + alpha^2*Vi2(a,b))) ;
        end            
    end
    
    [c,p] = corrcoef(sn_corrCoefs(1:numel(Pss_predicted))',Pss_predicted(:)) ; 
    CorrActPred_sn(iter) = c(1,2) ;
    Mse_sn(iter) = mean((sn_corrCoefs(1:numel(Pss_predicted)) - Pss_predicted(:)').^2) ;
    
    [c,p] = corrcoef(sn_corrCoefs(1:numel(Pss_predictedNoCross))',Pss_predictedNoCross(:)) ; 
    CorrActPredNoCross_sn(iter) = c(1,2) ;
    MseNoCross_sn(iter) = mean((sn_corrCoefs(1:numel(Pss_predictedNoCross)) - Pss_predictedNoCross(:)').^2) ;
    
    [c,p] = corrcoef(sT_corrCoefs(1:numel(Pss_predicted))',Pss_predicted(:)) ; 
    CorrActPred_st(iter) = c(1,2) ;
    Mse_st(iter) = mean((sT_corrCoefs(1:numel(Pss_predicted)) - Pss_predicted(:)').^2) ;
    
    [c,p] = corrcoef(sT_corrCoefs(1:numel(Pss_predictedNoCross))',Pss_predictedNoCross(:)) ; 
    CorrActPredNoCross_st(iter) = c(1,2) ;
    MseNoCross_st(iter) = mean((sT_corrCoefs(1:numel(Pss_predictedNoCross)) - Pss_predictedNoCross(:)').^2) ;
    
end


% preedict using the best alpha
% i1 = find(Mse_sn==min(Mse_sn)) ;
% alpha_sn = alpha_range(i1) ;
% 
% i1 = find(Mse_st==min(Mse_st)) ;
% alpha_st = alpha_range(i1) ;


% predict using alpha from Voltage_mean
alpha = abs((Voltage_mean-InhRev)/(Voltage_mean-0)) ; % alpha from voltage mean
alpha = round(alpha*10/(alpha_range(2)-alpha_range(1)))/(10/(alpha_range(2)-alpha_range(1))) ; % round to the nearest value within alpha range
errori = find(alpha_range==alpha) ; % the index of the error for the alpha taken from Voltage_mean  

alpha_sn = alpha ;
alpha_st = alpha ;

for a=1:numBars; % for every bar
    for b=1:numConditions ; % for every set pairings   
        numerator = Ve1e2(a,b) + alpha_sn^2*Vi1i2(a,b) - alpha_sn*Ve1i2(a,b) - alpha_sn*Vi1e2(a,b) ;
        denominator = sqrt((Ve1(a,b) + alpha_sn^2*Vi1(a,b) - 2*alpha_sn*Ve1i1(a,b))*(Ve2(a,b) + alpha_sn^2*Vi2(a,b) - 2*alpha_sn*Ve2i2(a,b))) ;
        Pss_predicted_sn(a,b) = numerator/denominator ;
        
        numerator = Ve1e2(a,b) + alpha_st^2*Vi1i2(a,b) - alpha_st*Ve1i2(a,b) - alpha_st*Vi1e2(a,b) ;
        denominator = sqrt((Ve1(a,b) + alpha_st^2*Vi1(a,b) - 2*alpha_st*Ve1i1(a,b))*(Ve2(a,b) + alpha_st^2*Vi2(a,b) - 2*alpha_st*Ve2i2(a,b))) ;
        Pss_predicted_st(a,b) = numerator/denominator ;
    end            
end



% % figures
% 
% figure % estimated alpha values
% plot(alphaHistX,alpha_snHist,'o-')
% hold on
% plot(alphaHistX,alpha_stHist,'+-')
% xlabel('alpha')
% ylabel('number of observations')
% 
% 
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
plot(sn_corrCoefs(1:numel(Pss_predicted)),Pss_predicted_sn(:),'o')
hold on
plot(sT_corrCoefs(1:numel(Pss_predicted)),Pss_predicted_st(:),'+')
xlabel('actual')
ylabel('predicted')
legend('sn','st')

% FOR IGOR
%{
identifier = ['Ve1e2',id,'cell',num2str(A)] ; % Ve1e2 
ForIgor.(identifier) = Ve1e2 ; 

identifier = ['Vi1i2',id,'cell',num2str(A)] ; % Vi1i2 
ForIgor.(identifier) = Vi1i2 ; 

identifier = ['Ve1i2',id,'cell',num2str(A)] ; % Ve1i2 
ForIgor.(identifier) = Ve1i2 ; 

identifier = ['Vi1e2',id,'cell',num2str(A)] ; % Vi1e2 
ForIgor.(identifier) = Vi1e2 ;

identifier = ['Ve1i1',id,'cell',num2str(A)] ; % Ve1i1
ForIgor.(identifier) = Ve1i1 ; 

identifier = ['Ve2i2',id,'cell',num2str(A)] ; % Ve2i2
ForIgor.(identifier) = Ve2i2 ; 

identifier = ['Ve1',id,'cell',num2str(A)] ; % Ve1
ForIgor.(identifier) = Ve1 ; 

identifier = ['Vi1',id,'cell',num2str(A)] ; % Vi1
ForIgor.(identifier) = Vi1 ; 

identifier = ['Ve2',id,'cell',num2str(A)] ; % Ve2
ForIgor.(identifier) = Ve1 ; 

identifier = ['Vi2',id,'cell',num2str(A)] ; % Vi2
ForIgor.(identifier) = Vi2 ; 

%}

% identifier = ['snCorrsActual',id,'cell',num2str(A)] ;
% ForIgor.(identifier) = sn_corrCoefs(1:numel(Pss_predicted)) ;
% 
% identifier = ['snCorrsAnalytic',id,'cell',num2str(A)] ;
% ForIgor.(identifier) = Pss_predicted_sn(:) ;
%
% identifier = ['sTCorrsActual',id,'cell',num2str(A)] ;
% ForIgor.(identifier) = sT_corrCoefs(1:numel(Pss_predicted)) ;
% 
% identifier = ['sTCorrsAnalytic',id,'cell',num2str(A)] ;
% ForIgor.(identifier) = Pss_predicted_st(:) ;

for a = 1:numBars ;
    identifier = ['linAnalyticsnCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = Pss_predicted_sn(a,:) ;     
end
    
    
%}
% make room in memmory
clearvars -except ForIgor Parameters id A excg* inhg* sT_corrCoefs sn_corrCoefs num* prePnts barPnts minusPairwiseShift SI time*


%% noise correlation coefs in presented conductances  

% g correlations
%{
for b=1:numBars ; % preallocate variables to improve memory usage and speed
    cc1_plusConverge{b} = nans(numTrials-2,barPnts*2-1) ;
    cc2_plusConverge{b} = nans(numTrials-2,barPnts*2-1) ;
    cc1_minusConverge{b} = nans(numTrials-2,barPnts*2-1) ;
    cc2_minusConverge{b} = nans(numTrials-2,barPnts*2-1) ;
    ccExc_plusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    ccInh_plusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    ccExc_minusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    ccInh_minusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    ccExc_plusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    ccInh_plusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    ccExc_minusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    ccInh_minusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;  
    ccExcInh_plusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    ccInhExc_plusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;     
    ccExcInh_minusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    ccInhExc_minusConverg_plusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    ccExcInh_plusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    ccInhExc_plusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    ccExcInh_minusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;
    ccInhExc_minusConverg_minusPairwise{b} = nans(numTrials-2,barPnts*2-1) ;  
end
cc1_plusConverge_peak = nans(numBars,numTrials-2) ;
cc2_plusConverge_peak = nans(numBars,numTrials-2) ;        
cc1_minusConverge_peak = nans(numBars,numTrials-2) ;
cc2_minusConverge_peak = nans(numBars,numTrials-2) ;
ccExc_plusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ;
ccInh_plusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ;       
ccExc_minusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ;
ccInh_minusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ;        
ccExc_plusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ;
ccInh_plusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ;       
ccExc_minusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ;
ccInh_minusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ;      
ccExcInh_plusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ; 
ccInhExc_plusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ; 
ccExcInh_minusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ; 
ccInhExc_minusConverg_plusPairwise_peak = nans(numBars,numTrials-2) ;         
ccExcInh_plusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ; 
ccInhExc_plusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ;        
ccExcInh_minusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ; 
ccInhExc_minusConverg_minusPairwise_peak = nans(numBars,numTrials-2) ;      
ccSpikePred_plusConverg_plusPairwise = nans(numBars,numTrials-2) ; 
ccSpikePred_minusConverg_plusPairwise = nans(numBars,numTrials-2) ; 
ccSpikePred_plusConverg_minusPairwise = nans(numBars,numTrials-2) ; 
ccSpikePred_minusConverg_minusPairwise = nans(numBars,numTrials-2) ; 
cc1_plusConverge_mean = nans(numBars,barPnts*2-1) ; 
cc2_plusConverge_mean = nans(numBars,barPnts*2-1) ;    
cc1_minusConverge_mean = nans(numBars,barPnts*2-1) ; 
cc2_minusConverge_mean = nans(numBars,barPnts*2-1) ;    
ccExc_plusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
ccInh_plusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
ccExc_minusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
ccInh_minusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ;    
ccExc_plusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ; 
ccInh_plusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ; 
ccExc_minusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ; 
ccInh_minusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ;   
ccExcInh_plusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
ccInhExc_plusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
ccExcInh_minusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
ccInhExc_minusConverg_plusPairwise_mean = nans(numBars,barPnts*2-1) ; 
ccExcInh_plusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ; 
ccInhExc_plusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ; 
ccExcInh_minusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ; 
ccInhExc_minusConverg_minusPairwise_mean = nans(numBars,barPnts*2-1) ;  


for b=1:numBars ;
    excg2_sim_bar_res_circshift{b} = circshift(excg2_sim_bar_res{b},[minusPairwiseShift,0]) ;
    inhg2_sim_bar_res_circshift{b} = circshift(inhg2_sim_bar_res{b},[minusPairwiseShift,0]) ;
    excg2_shuff_bar_res_circshift{b} = circshift(excg2_shuff_bar_res{b},[minusPairwiseShift,0]) ;
    inhg2_shuff_bar_res_circshift{b} = circshift(inhg2_shuff_bar_res{b},[minusPairwiseShift,0]) ;
    
    for a = 1:numTrials-2 ; % for each residual
        
        % converging correlations coef
        cc1_plusConverge{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),inhg1_sim_bar_res{b}(a,:),'coef') ;
        cc2_plusConverge{b}(a,:) = xcorr(excg2_sim_bar_res{b}(a,:),inhg2_sim_bar_res{b}(a,:),'coef') ;
        
        cc1_minusConverge{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),inhg1_shuff_bar_res{b}(a,:),'coef') ;
        cc2_minusConverge{b}(a,:) = xcorr(excg2_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res{b}(a,:),'coef') ;
        
    
        % pairwise ee,ii correlations coefs
        ccExc_plusConverg_plusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),excg2_sim_bar_res{b}(a,:),'coef') ; % + converging + pairwise
        ccInh_plusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res{b}(a,:),'coef') ;
        
        ccExc_minusConverg_plusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res{b}(a,:),'coef') ; % -converging + pairwise
        ccInh_minusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res{b}(a,:),'coef') ;
        
        ccExc_plusConverg_minusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),excg2_sim_bar_res_circshift{b}(a,:),'coef') ; % + converging - pairwise
        ccInh_plusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res_circshift{b}(a,:),'coef') ;
        
        ccExc_minusConverg_minusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res_circshift{b}(a,:),'coef') ; % -converging - pairwise
        ccInh_minusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res_circshift{b}(a,:),'coef') ;  

        % pairwise ei,ie correlations coefs
        ccExcInh_plusConverg_plusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res{b}(a,:),'coef') ; % + converging + pairwise
        ccInhExc_plusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),excg2_sim_bar_res{b}(a,:),'coef') ;
        
        ccExcInh_minusConverg_plusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res{b}(a,:),'coef') ; % -converging + pairwise
        ccInhExc_minusConverg_plusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res{b}(a,:),'coef') ;
        
        ccExcInh_plusConverg_minusPairwise{b}(a,:) = xcorr(excg1_sim_bar_res{b}(a,:),inhg2_sim_bar_res_circshift{b}(a,:),'coef') ; % + converging - pairwise
        ccInhExc_plusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_sim_bar_res{b}(a,:),excg2_sim_bar_res_circshift{b}(a,:),'coef') ;
        
        ccExcInh_minusConverg_minusPairwise{b}(a,:) = xcorr(excg1_shuff_bar_res{b}(a,:),inhg2_shuff_bar_res_circshift{b}(a,:),'coef') ; % -converging - pairwise
        ccInhExc_minusConverg_minusPairwise{b}(a,:) = xcorr(inhg1_shuff_bar_res{b}(a,:),excg2_shuff_bar_res_circshift{b}(a,:),'coef') ;  
  
        
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

% PCA (the goal here is to solve the output correlation as a function of the wieghted input
% correlations: ie. a5*Pss = a1*Pee + a2*Pii + a3*Pei + a4*Pie ; the vector
% a is the princple component)

% every bar has a spike number correlation and mean g correlation
Pee = [ccExc_plusConverg_plusPairwise_mean_peak(:),ccExc_minusConverg_plusPairwise_mean_peak(:),ccExc_plusConverg_minusPairwise_mean_peak(:),ccExc_minusConverg_minusPairwise_mean_peak(:)] ;
Pii = [ccInh_plusConverg_plusPairwise_mean_peak(:),ccInh_minusConverg_plusPairwise_mean_peak(:),ccInh_plusConverg_minusPairwise_mean_peak(:),ccInh_minusConverg_minusPairwise_mean_peak(:)] ;
Pei = [ccExcInh_plusConverg_plusPairwise_mean_peak(:),ccExcInh_minusConverg_plusPairwise_mean_peak(:),ccExcInh_plusConverg_minusPairwise_mean_peak(:),ccExcInh_minusConverg_minusPairwise_mean_peak(:)] ;
Pie = [ccInhExc_plusConverg_plusPairwise_mean_peak(:),ccInhExc_minusConverg_plusPairwise_mean_peak(:),ccInhExc_plusConverg_minusPairwise_mean_peak(:),ccInhExc_minusConverg_minusPairwise_mean_peak(:)] ;

M = [Pee(:),Pii(:),Pei(:),Pie(:),sn_corrCoefs(:)] ;

i1 = find(isnan(M)==1) ;
M(i1) = 0 ;

Mcov = cov(M) ;

[EigVec,EigVal] = eig(Mcov) ;

Variances_sn = 100*diag(EigVal)/sum(diag(EigVal)) ;

i1 = find(Variances_sn == max(Variances_sn)) ;
PrinEigVec_sn = EigVec(:,i1) ;

Pss_estimate_sn = (M(:,1:4)*PrinEigVec_sn(1:4))/PrinEigVec_sn(5) ;
Frac_variancePss_sn = 1- mean((sn_corrCoefs(:) - Pss_estimate_sn(:)).^2)/mean((sn_corrCoefs(:) - mean(sn_corrCoefs(:))).^2) ; % fraction of variance of spike output predicted by coefficients
temp = corrcoef(sn_corrCoefs(:), Pss_estimate_sn(:)) ;
CorrPssEstPssTrue_sn = temp(2,1) ;

M = [Pee(:),Pii(:),Pei(:),Pie(:),sT_corrCoefs(:)] ;

i1 = find(isnan(M)==1) ;
M(i1) = 0 ;

Mcov = cov(M) ;

[EigVec,EigVal] = eig(Mcov) ;

Variances_st = 100*diag(EigVal)/sum(diag(EigVal)) ;

i1 = find(Variances_sn == max(Variances_sn)) ;
PrinEigVec_st = EigVec(:,i1) ;

Pss_estimate_st = (M(:,1:4)*PrinEigVec_sn(1:4))/PrinEigVec_sn(5) ;
Frac_variancePss_st = 1- mean((sn_corrCoefs(:) - Pss_estimate_st(:)).^2)/mean((sn_corrCoefs(:) - mean(sn_corrCoefs(:))).^2) ; % fraction of variance of spike output predicted by coefficients
temp = corrcoef(sn_corrCoefs(:), Pss_estimate_st(:)) ;
CorrPssEstPssTrue_st = temp(2,1) ;


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

% input correlations
for a = 1:numBars ;
    identifier = ['Pe1e2',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = Pee(a,:);

    identifier = ['Pi1i2',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = Pii(a,:);
    
    identifier = ['Pe1i2',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = Pei(a,:);
    
    identifier = ['Pi1e2',id,'Bar',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = Pie(a,:);
   
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