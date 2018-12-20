function ForIgor = DCdsPairRepeatAnalyzer(Input,Parameters,id,A) ;

% this script is a modifiead version of DCdsPairAnalyzer.m and is designed
% to assess the extent of correlation between a repeated sequence a injected condtances.  
% Input data should be organized sets 1,2,3,4,1,2,3,4  

% JC 3/23/11

epochs = str2num(Input(A).(id)) ; %sim and shuffled g

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);

for a = 1:length(epochs) ;
    [data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ; 
    [excg(a,:), inhg(a,:), error(a,:)] = ITCReadEpochStmGClamp(epochs(a), 0, fp);
end

% get sampling interval 
[SI, error] = ITCGetSamplingInterval(epochs(1), fp);    
SI = SI * 1e-6; % Sampling interval in sec
if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end
time = [SI:SI:SI*length(data)] ;

% spike detection
SpikePnts = SpikeDetection_WC(data,-30,10000) ; 

% spike trains
numTrials = length(SpikePnts)/2 ; % number of mimiced trials for each set of repeated g

spikeTrain1 = zeros(numTrials,length(data)) ; 
spikeTrain2 = zeros(numTrials,length(data)) ; 

for a = 1:numTrials ;
    spikeTrain1(a,SpikePnts{a})=1 ; % first time g are presented
    spikeTrain2(a,SpikePnts{numTrials+a})=1 ; % repetition of g 
end

% spike trains divided by bar
numBars = 3 ; % number of bars presented while recording conductances
prePnts = 2000 ; % pnts to avoid as g starts to be injected
barPnts = floor((length(data)-prePnts)/numBars) ; 

for a=1:numBars ;
    spikeTrain1_bar{a} = spikeTrain1(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
    spikeTrain2_bar{a} = spikeTrain2(:,prePnts+barPnts*a-barPnts+1:prePnts+barPnts*a) ;
end

% spike number for each bar
for a=1:numBars ;
    sn1(:,a) = sum(spikeTrain1_bar{a},2) ;
    sn2(:,a) = sum(spikeTrain2_bar{a},2) ;

    sn1_mean(a) = mean(sn1(:,a)) ;
    sn2_mean(a) = mean(sn2(:,a)) ;

end

% smoothed spike trains
smthFctr = .050 ; % sec
smthFctr = round(smthFctr/SI) ;
for a=1:numBars ;
    for b=1:numTrials ;
        spikeTrain1_smth_bar{a}(b,:) = smooth(spikeTrain1_bar{a}(b,:),smthFctr) ;
        spikeTrain2_smth_bar{a}(b,:) = smooth(spikeTrain2_bar{a}(b,:),smthFctr) ;
    end
end

% spike number residuals (local correction)
for a = 2:numTrials-1 ;
    sn1_res(a-1,:) = sn1(a,:) - (sn1(a-1,:)+sn1(a+1,:))/2 ;
    sn2_res(a-1,:) = sn2(a,:) - (sn2(a-1,:)+sn2(a+1,:))/2 ;
end

% smoothed spike train residuals
for a=1:numBars ;
    for b=2:numTrials-1  ;
        spikeTrain1_smth_bar_res{a}(b-1,:) = spikeTrain1_smth_bar{a}(b,:) - (spikeTrain1_smth_bar{a}(b-1,:) + spikeTrain1_smth_bar{a}(b+1,:))./2 ;
        spikeTrain2_smth_bar_res{a}(b-1,:) = spikeTrain2_smth_bar{a}(b,:) - (spikeTrain2_smth_bar{a}(b-1,:) + spikeTrain2_smth_bar{a}(b+1,:))./2 ;  
    end
end

% sn unbiased variance of spike number from residuals
sn1_var = sum(sn1_res.^2,1)./(size(sn1_res,1)-1) ;
sn2_var = sum(sn2_res.^2,1)./(size(sn2_res,1)-1) ;

% spike number correlation (linear) for each bar
lincoefsigMin = .05 ;
for a=1:numBars ;
    [Corr,p] =  corrcoef(sn1_res(:,a),sn2_res(:,a)) ; % correlations 
    sn_repeated_lincoef(a) = Corr(1,2) ;
    sn_repeated_lincoefsig(a) = p(1,2) ;
    if sn_repeated_lincoefsig(a)>lincoefsigMin ;
        sn_repeated_lincoef(a) = 0 ;
    end
end

% smoothed spike train correlations
for a=1:numBars ;
   
    for b=1:numTrials-2  ;
        sT_repeated_cc{a}(b,:) = xcov(spikeTrain1_smth_bar_res{a}(b,:),spikeTrain2_smth_bar_res{a}(b,:),'coef') ;

        % peaks of cc
        sT_repeated_ccPeak(a,b) = CCpeakFinder(sT_repeated_cc{a}(b,:)) ;

    end
    
    sT_repeated_cc_mean(a,:) = nanmean(sT_repeated_cc{a}) ; % trials where residual = 0 have no cc (no residual = no noise = no ability to correlate noise)
    
    sT_repeated_ccPeak_mean(a) = nanmean(sT_repeated_ccPeak(a,:)) ; % average of peaks not peak of average 

end
    
% conductances
for a = 1:numTrials ;
    excg1(a,:) = excg(a,:) ; 
    excg2(a,:) = excg(numTrials+a,:) ;  

    inhg1(a,:) = inhg(a,:) ; 
    inhg2(a,:) = inhg(numTrials+a,:) ; 
end

if excg1~=excg2 & inhg1~=inhg2 ;
    disp('not repeated g')
end

time_cc = SI*([1:length(sT_repeated_cc_mean)]-(length(sT_repeated_cc_mean)+1)/2) ;

% figures
figure % spike detection check
for a=1:numTrials ; % on each trial
    subplot(2,1,1)
    plot(time,data(a,:))
    hold on
    plot(time,spikeTrain1(a,:)*-30,'r')
    hold off
    
    text(0,.9,num2str(epochs(a)),'units','norm')
    
    subplot(2,1,2)
    plot(time,data(a+numTrials,:))
    hold on
    plot(time,spikeTrain2(a,:)*-30,'r')
    hold off    
    
    text(0,.9,num2str(epochs(a+numTrials)),'units','norm')
    pause
end

figure % spike number noise correlations
for a = 1:numBars ;
    subplot(1,numBars,a)
    plot(sn1_res(:,a),sn2_res(:,a),'k*')  
end

figure % spike train noise correlations
for a = 1:numBars ;
subplot(1,numBars,a)
    plot(time_cc,sT_repeated_cc_mean(a,:),'k') 


    xlim([-1 1])
    ylim([-.5 1])
end

% FOR IGOR

for a = 1:numBars ;
    identifier = ['snCorrsBar','id',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sn_repeated_lincoef(a) ;
    
    identifier = ['stCorrsBar','id',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = sT_repeated_ccPeak_mean(a) ;    
end



