function ForIgor = pairCADS(Input,Parameters,id,A) ; 
% this function will analyze a pair of ds cells in responses to a moving
% bar

% 8/30/10 added break points to calculate local mean and made so that epochs could be different lengths 
% 2/2/11 added MonitorAnalysis and poisson export stuff

% get data
[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
numTrials = length(epochs) 

for a = 1:numTrials ; % for each spike epoch
    [dataCA1{a}, error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data
    [dataCA2{a}, error] = ITCReadEpoch(epochs(a), 1, fp) ;    % get data
    
    numPnts(a) = length(dataCA1{a}) ;
    
    [SI_CA(a), error] = ITCGetSamplingInterval(epochs(a), fp); % get sampling interval
    SI_CA(a) = SI_CA(a) * 1e-6; % Sampling interval in sec
    
    [epochstarttime(a), error] = ITCGetEpochTime(epochs(a),fp) ; % seconds since start of cell file
end

if Input(A).ITC18flag == 1 ;
    SI_CA = SI_CA*1.25 ;
end

SpatialStimParams = hdf5load(['~/Data/mouse/',Input(A).cellname,'_spatial.h5']) ; % load spatial stim params

Monitor = nans(numTrials,max(numPnts)) ;
for a = 1:numTrials ;
    [Monitor(a,1:numPnts(a)), error] = ITCReadEpoch(epochs(a), 5, fp) ;    % get monitor data
end
[framerate_mean,framerate_range,startAmp] = MonitorAnalysis(Monitor, length(epochs), SI_CA(1)) ; % get frame rate
frameRate = mean(framerate_mean) ;
% frameRate = 60 ;

% get spatial stim parameters
BarAngles_CA_unique = [] ;

for a = 1:numTrials ;
    time{a} = [1:length(dataCA1{a})]*SI_CA(a) ; % time for epoch
    
    % get spatial stim parameters
    StrucString = ['params_epoch',num2str(epochs(a))] ; 
    Struct = SpatialStimParams.(StrucString) ;

    BarAngles_CA{a} = Struct.BarAngle ;
    BarWidth_CA(a) = Struct.BarWidth ;
    BarSpeed_CA(a) = Struct.BarSpeed ;
    
    prePnts_CA(a) = floor(Struct.spatial_prepts/(frameRate*SI_CA(a))) ;
    postPnts_CA(a) = floor(Struct.spatial_postpts/(frameRate*SI_CA(a))) ;
    groupPnts_CA(a) = floor((length(dataCA1{a})-prePnts_CA(a)-postPnts_CA(a))/length(BarAngles_CA{a})) ;
    
    NumBars_CA(a) = length(BarAngles_CA{a}) ; % number of bars shown per trial
    BarPnts_CA(a) = floor((length(dataCA1{a})-prePnts_CA(a)-postPnts_CA(a))/NumBars_CA(a)) ; % number of points the bar is presented + interbar points
    OnPnts_CA(a) = floor((BarWidth_CA(a)/BarSpeed_CA(a))/(frameRate*SI_CA(1))) ; % number of points before end of bar appears triggering off response

    timeBar{a} = [1:BarPnts_CA(a)]*SI_CA(a) ; % time of each bar
    
    BarAngles_CA_unique = unique([BarAngles_CA{a};BarAngles_CA_unique]) ; % unique bar angles
end
numUniqueBars = length(BarAngles_CA_unique) ;

% prep NaN matricies (each row is a trial, each column is a particular bar angle)
SN1 = nans(numTrials,numUniqueBars) ;
SN2 = SN1 ;
ST1 = SN1 ;
ST2 = SN1 ;
SN1_res = SN1 ;
SN2_res = SN1 ;
ST1_res = SN1 ;
ST2_res = SN1 ;

% get spike stats
spikeThreshold1 = 15;
spikeThreshold2 = 15 ;

for a = 1:numTrials ;
    % spike detections
    SpikePnts1{a} = SpikeDetection(dataCA1{a},spikeThreshold1,1/SI_CA(a)) ;
    SpikePnts2{a} = SpikeDetection(dataCA2{a},spikeThreshold2,1/SI_CA(a)) ;

    SpikeTrain1{a} = zeros(size(dataCA1{a})) ;
    SpikeTrain2{a} = zeros(size(dataCA2{a})) ;

    SpikeTrain1{a}(SpikePnts1{a}{1}) = 1 ;
    SpikeTrain2{a}(SpikePnts2{a}{1}) = 1 ;

    % spike number for each bar presented on this trial
    for  b = 1:NumBars_CA(a) ; % for each bar
        unique_i = find(BarAngles_CA_unique==BarAngles_CA{a}(b)) ;
        
        SpikeTrainBar1{a}{unique_i} = SpikeTrain1{a}(prePnts_CA(a)+BarPnts_CA(a)*(b-1):prePnts_CA(a)+BarPnts_CA(a)*b) ; % spike train for each bar
        SpikeTrainBar2{a}{unique_i} = SpikeTrain2{a}(prePnts_CA(a)+BarPnts_CA(a)*(b-1):prePnts_CA(a)+BarPnts_CA(a)*b) ;
        
        timeCheck{a}(b,:) = time{a}(prePnts_CA(a)+BarPnts_CA(a)*(b-1):prePnts_CA(a)+BarPnts_CA(a)*b) ;
        
        SN1(a,unique_i) = sum(SpikeTrainBar1{a}{unique_i}) ; % spike number
        SN2(a,unique_i) = sum(SpikeTrainBar2{a}{unique_i}) ;
        
        % time of first spike
        if SN1(a,unique_i)>0 ;
            ST1(a,unique_i) = timeBar{a}(find(SpikeTrainBar1{a}{unique_i}==1,1)) ;
        end
        
        if SN2(a,unique_i)>0 ;
            ST2(a,unique_i) = timeBar{a}(find(SpikeTrainBar2{a}{unique_i}==1,1)) ;
        end
    end
end

% difference in time of first spikes
STdiff = ST1-ST2 ; 

% get resduals where appropriate (where sequential preceding and following epochs exist)
breakPntsSecs = 5 ; % acceptable difference between trial start times where epochs are still considered sequential
epochsSecs = diff(epochstarttime) ;

for a = 2:numTrials-1 ;
    if abs(epochsSecs(a)-epochsSecs(a-1))<breakPntsSecs & NumBars_CA(a)==NumBars_CA(a-1) & NumBars_CA(a)==NumBars_CA(a+1) ; % if preceding and following epochs are sequential then calculate a residual
        SN1_res(a,:) = SN1(a,:) - ((SN1(a-1,:)+SN1(a+1,:))./2) ;
        SN2_res(a,:) = SN2(a,:) - ((SN2(a-1,:)+SN2(a+1,:))./2) ;
        
        ST1_res(a,:) = ST1(a,:) - ((ST1(a-1,:)+ST1(a+1,:))./2) ;
        ST2_res(a,:) = ST2(a,:) - ((ST2(a-1,:)+ST2(a+1,:))./2) ;
    end
end

% get correlation of residuals
resNumTrialsMin = 20 ; % the minimum number of trials for which we will still try to calculate a correlation coef
lincoefsigMin = .05 ; % below this sig value corr values will be recorded as zero
numStd = 2 ; % number of std outside mean from which trials are calculated for corr coef

for a=1:numUniqueBars ;
    SN_i = find(~isnan(SN1_res(:,a))) ; % use trials where a residual was calculated
    
    SN1_var(a) = sum(SN1_res(SN_i,a).^2)/length(SN_i) ; % calculate biased variance from residuals which are locally corrected for drift
    SN2_var(a) = sum(SN2_res(SN_i,a).^2)/length(SN_i) ;
    
    if length(SN_i)>resNumTrialsMin ; % if there are enough then calculate a corr coef
        SN1_res_group{a} = SN1_res(SN_i,a) ; % residual values
        SN2_res_group{a} = SN2_res(SN_i,a) ;
        
        SN1_res_mean = mean(SN1_res_group{a}) ; % mean of residuals
        SN2_res_mean = mean(SN2_res_group{a}) ; 
        
        SN1_res_std = std(SN1_res_group{a}) ; % std of residuals
        SN2_res_std = std(SN2_res_group{a}) ; 
        
        SN1_res_group_i = find(SN1_res_group{a}>SN1_res_mean-numStd*SN1_res_std & SN1_res_group{a}<SN1_res_mean+numStd*SN1_res_std) ; % those trials that are within x std of mean
        SN2_res_group_i = find(SN2_res_group{a}>SN2_res_mean-numStd*SN2_res_std & SN2_res_group{a}<SN2_res_mean+numStd*SN2_res_std) ;
        SN_res_group_i{a} = intersect(SN1_res_group_i,SN2_res_group_i) ; % only use trials where both residuals are with 2 std of mean
        SN_resNumTrials(a) = length(SN_res_group_i{a}) ; % number of trials used to calculate correlation coef
        
       
%         SN_res_lincoef(a) = sum(SN1_res_group(SN_res_group_i).*SN2_res_group(SN_res_group_i))/((SN_resNumTrials(a)-1)*(std(SN1_res_group(SN_res_group_i))*std(SN2_res_group(SN_res_group_i)))) ; % this gives a slightly diff number than matlabs function?
     
        [Corr,p] =  corrcoef(SN1_res_group{a}(SN_res_group_i{a}),SN2_res_group{a}(SN_res_group_i{a})) ;  % linear correlation (pearsons)
        SN_res_lincoef(a) = Corr(1,2) ;
        SN_res_lincoefsig(a) = p(1,2) ;
        if SN_res_lincoefsig(a)>lincoefsigMin ;
            SN_res_lincoef(a) = 0 ;
        end
    else
        SN1_res_group{a} = nan ;
        SN2_res_group{a} = nan ;
        SN_res_group_i{a} = nan ;
        SN_resNumTrials(a) = nan ;
        SN_res_lincoef(a) = nan ;
        SN_res_lincoefsig(a) = nan ;  
    end
    
    ST_i = intersect(find(~isnan(ST1_res(:,a))),find(~isnan(ST2_res(:,a)))) ; % use trials only where both cells spiked
    if length(ST_i)>resNumTrialsMin ;  % if there are enough then calculate a corr coef
        ST1_res_group{a} = ST1_res(ST_i,a) ; % residual values
        ST2_res_group{a} = ST2_res(ST_i,a) ;
        
        ST1_res_mean = mean(ST1_res_group{a}) ; % mean of residuals
        ST2_res_mean = mean(ST2_res_group{a}) ; 
        
        ST1_res_std = std(ST1_res_group{a}) ; % std of residuals
        ST2_res_std = std(ST2_res_group{a}) ; 
        
        ST1_res_group_i = find(ST1_res_group{a}>ST1_res_mean-numStd*ST1_res_std & ST1_res_group{a}<ST1_res_mean+numStd*ST1_res_std) ; % those trials that are within 2 std of mean
        ST2_res_group_i = find(ST2_res_group{a}>ST2_res_mean-numStd*ST2_res_std & ST2_res_group{a}<ST2_res_mean+numStd*ST2_res_std) ;
        ST_res_group_i{a} = intersect(ST1_res_group_i,ST2_res_group_i) ; % only use trials where both residuals are with 2 std of mean
        ST_resNumTrials(a) = length(ST_res_group_i{a}) ; % number of trials used to calculate correlation coef
        
        [Corr,p] =  corrcoef(ST1_res_group{a}(ST_res_group_i{a}),ST2_res_group{a}(ST_res_group_i{a})) ; % linear correlation (pearsons)
        ST_res_lincoef(a) = Corr(1,2) ;
        ST_res_lincoefsig(a) = p(1,2) ;
        if ST_res_lincoefsig(a)>lincoefsigMin ;
            ST_res_lincoef(a) = 0 ;
        end
    else 
        ST1_res_group{a} = nan ;
        ST2_res_group{a} = nan ;
        ST_res_group_i{a} = nan ;
        ST_resNumTrials(a) = nan ;
        ST_res_lincoef(a) = nan ;
        ST_res_lincoefsig(a) = nan ;
    end
    
end
 
% get spike distrubtions
for a=1:numUniqueBars ;
  
    [SN1_dist{a},SN1_distX{a}] = hist(SN1(:,a),[min(SN1(:,a)):max(SN1(:,a))]) ; 
    [SN2_dist{a},SN2_distX{a}] = hist(SN2(:,a),[min(SN2(:,a)):max(SN2(:,a))]) ;
    
    STbin = .01 ; % sec
    [ST1_dist{a},ST1_distX{a}] = hist(ST1(:,a),[min(ST1(:,a))-STbin:STbin:max(ST1(:,a))]) ; 
    [ST2_dist{a},ST2_distX{a}] = hist(ST2(:,a),[min(ST2(:,a))-STbin:STbin:max(ST2(:,a))]) ; 
    
    STdiffbin = .01 ; % sec
    [STdiff_dist{a},STdiff_distX{a}] = hist(STdiff(:,a),[min(STdiff(:,a))-STdiffbin:STdiffbin:max(STdiff(:,a))]) ; 

    SN1_std(a) = nanstd(SN1(:,a)) ; % FIX THESE STD THEY SHOULD BE CALC FROM LOCAL CORRECTED RESID
    SN2_std(a) = nanstd(SN2(:,a)) ;

    ST1_std(a) = nanstd(ST1(:,a)) ;
    ST2_std(a) = nanstd(ST2(:,a)) ;  
    STdiff_std(a) = nanstd(STdiff(:,a)) ;
    
    SN1_mean(a) = nanmean(SN1(:,a)) ;
    SN2_mean(a) = nanmean(SN2(:,a)) ;

    ST1_mean(a) = nanmean(ST1(:,a)) ;
    ST2_mean(a) = nanmean(ST2(:,a)) ; 
    STdiff_mean(a) = nanmean(STdiff(:,a)) ;

end
    
% raster prep
for  a = 1:numTrials ; % for each trial
    for b = 1:NumBars_CA(a) ;
        unique_i = find(BarAngles_CA_unique==BarAngles_CA{a}(b)) ;
        
        raster1{unique_i}{a} = timeBar{a}(SpikeTrainBar1{a}{unique_i}==1) ;
        raster2{unique_i}{a} = timeBar{a}(SpikeTrainBar2{a}{unique_i}==1) ;
    end
end

% vector tuning 
for a=1:numUniqueBars ; %for each bar
    
    vector_magX_1(a) = SN1_mean(a)*cosd(BarAngles_CA_unique(a)) ; % magnitude of individual vector components
    vector_magY_1(a) = SN1_mean(a)*sind(BarAngles_CA_unique(a)) ;
       
    vector_magX_2(a) = SN2_mean(a)*cosd(BarAngles_CA_unique(a)) ;
    vector_magY_2(a) = SN2_mean(a)*sind(BarAngles_CA_unique(a)) ;
           
end

vectorMagnitude1 = sqrt(sum(vector_magX_1)^2+sum(vector_magY_1)^2) ; % magnitude of summed vector
vectorMagnitude2 = sqrt(sum(vector_magX_2)^2+sum(vector_magY_2)^2) ;

vectorAngle1 = atand(sum(vector_magX_1)/sum(vector_magY_1)) ; % angle of summed vector
vectorAngle2 = atand(sum(vector_magX_2)/sum(vector_magY_2)) ;       

if sum(vector_magX_1)<0 & sum(vector_magY_1)>0 ; % correct for atand ingorning quadrant
    vectorAngle1 = vectorAngle1 + 180 ;
elseif sum(vector_magX_1)<0 & sum(vector_magY_1)<0 ;
    vectorAngle1 = vectorAngle1 + 180 ;
elseif sum(vector_magX_1)>0 & sum(vector_magY_1)<0 ;
    vectorAngle1 = vectorAngle1 + 360 ; 
end

if sum(vector_magX_2)<0 & sum(vector_magY_2)>0 ;
    vectorAngle2 = vectorAngle2 + 180 ;
elseif sum(vector_magX_2)<0 & sum(vector_magY_2)<0 ;
    vectorAngle2 = vectorAngle2 + 180 ;
elseif sum(vector_magX_2)>0 & sum(vector_magY_2)<0 ;
    vectorAngle2 = vectorAngle2 + 360 ; 
end
    
vectorAngleDiff = abs(vectorAngle1-vectorAngle2) ;
if vectorAngleDiff>180
    vectorAngleDiff = 360-vectorAngleDiff ;
end

% clear variables for memmory
clearvars -except time Monitor...
    numTrials spikeThreshold* dataCA* SpikeTrain* epochs...
    numUniqueBars BarAngles_CA_unique SN* ST*...
    raster* vector* Input A id 
    
    
    

%figures

% figure
% for a = 1:numTrials ;
%     plot(time{a},Monitor)
%     hold on
% end
% 
% 
figure % spike detection check
for a=1:numTrials ; % on each trial
    subplot(2,1,1)
    plot(time{a},dataCA1{a})
    hold on
    plot(time{a},SpikeTrain1{a}*10+spikeThreshold1,'r')
    hold off
    
    subplot(2,1,2)
    plot(time{a},dataCA2{a})
    hold on
    plot(time{a},SpikeTrain2{a}*10+spikeThreshold2,'r')
    hold off
    
    text(0,.9,num2str(epochs(a)),'units','norm')
    pause
end
% 
% 
% figure % spike number stability
% subplot(1,2,1)
% for a = 1:numUniqueBars ; % for each bar
%     plot([1:numTrials],SN1(:,a),'b')
%     hold on
%     plot([1:numTrials],SN2(:,a),'r')
% end
% xlabel('trial')
% ylabel('spike number')
% 
% subplot(1,2,2)  % in entire trial
% plot(nansum(SN1,2),'b')
% hold on
% plot(nansum(SN2,2),'r')
% xlabel('trial')
% ylabel('total spike number')
% 
% 
% figure % rastors
% for a=1:numUniqueBars ;
%     subplot(1,numUniqueBars,a)
%     for b = 1:length(raster1{a}) ;
%         plot1ras(raster1{a}{b},b,.5,'b',.1) ;
%         hold on
%         plot1ras(raster2{a}{b},b+.5,.5,'r',.1) ;
%     end
% end
% 
% 
% figure % tuning curve comparison
% plot(BarAngles_CA_unique,SN1_mean,'b')
% hold on
% plot(BarAngles_CA_unique,SN2_mean,'r')
% for  a = 1:numUniqueBars ; % for each bar
%     plot(BarAngles_CA_unique(a),SN1(:,a),'bo')
%     plot(BarAngles_CA_unique(a),SN2(:,a),'r*')
% end
% legend('cell 1','cell 2')
% xlabel('bar angle')
% ylabel('number spikes')
% 
% figure % polar plots
% polar(0,max([vectorMagnitude1,vectorMagnitude2,SN1_mean,SN2_mean]))
% hold on
% 
% polar([vectorAngle1,vectorAngle1]*pi/180,[0,vectorMagnitude1],'b-')
% polar([vectorAngle2,vectorAngle2]*pi/180,[0,vectorMagnitude2],'r-')
% 
% polar([0:45:315]*pi/180,SN1_mean,'b*')
% polar([0:45:315]*pi/180,SN2_mean,'r*')
% disp(['vector angle difference = ',num2str(vectorAngleDiff)])
% 
% 
% figure % noise correlations?
% numsigdigits = 2 ;
% for  a = 1:numUniqueBars ; % for each bar
%     subplot(2,numUniqueBars,a)
%     plot(SN1_res(:,a),SN2_res(:,a),'k*')
%     hold on
%     if ~isnan(SN1_res_group{a})
%         plot(SN1_res_group{a}(SN_res_group_i{a}),SN2_res_group{a}(SN_res_group_i{a}),'y*') 
%     end
%     
%     %text(0,.9,['cc(sig)= ',num2str(SN_res_lincoef(a),numsigdigits),'(',num2str(SN_res_lincoefsig(a),numsigdigits),')'],'units','norm')
%     text(0,.9,['cc(sig)= ',num2str(SN_res_lincoef(a),numsigdigits)],'units','norm')
%     text(0,.8,['numtrials= ',num2str(SN_resNumTrials(a))],'units','norm')
%     
%     subplot(2,numUniqueBars,a+numUniqueBars)
%     plot(ST1_res(:,a),ST2_res(:,a),'k*')
%     hold on
%     if ~isnan(ST1_res_group{a})
%         plot(ST1_res_group{a}(ST_res_group_i{a}),ST2_res_group{a}(ST_res_group_i{a}),'y*') 
%     end
%     
%     %text(0,.9,['cc(sig)= ',num2str(ST_res_lincoef(a),numsigdigits),'(',num2str(ST_res_lincoefsig(a),numsigdigits),')'],'units','norm')
%     text(0,.9,['cc(sig)= ',num2str(ST_res_lincoef(a),numsigdigits)],'units','norm')
%     text(0,.8,['numtrials= ',num2str(ST_resNumTrials(a))],'units','norm')
% end
% 
% 
% figure % spike number and 1st spike time distributions cell 1
% for  a = 1:numUniqueBars ; % for each bar
%     if ~isempty(SN1_dist{a}) ;
%          subplot(3,numUniqueBars,a)
%          plot(SN1_distX{a},SN1_dist{a})
%          hold on
%          
%          subplot(3,numUniqueBars,a+numUniqueBars)
%          plot(ST1_distX{a},ST1_dist{a})
%          hold on
% 
%     end
%     if ~isempty(SN2_dist{a}) ;
%          subplot(3,numUniqueBars,a)
%          plot(SN2_distX{a},SN2_dist{a},'r')
%          
%          subplot(3,numUniqueBars,a+numUniqueBars)
%          plot(ST2_distX{a},ST2_dist{a},'r')
% 
%     end    
%     if ~isempty(STdiff_dist{a}) ;
%          subplot(3,numUniqueBars,a+2*numUniqueBars)
%          plot(STdiff_distX{a},STdiff_dist{a},'k')
% 
%     end    
% end

figure % mean spike number vs. var of spike number (is it poisson?)
unityVector = [min([SN1_mean,SN2_mean]),max([SN1_mean,SN2_mean])] ;
plot(SN1_mean,SN1_var,'*')
hold on
plot(SN2_mean,SN2_var,'r*')
plot(unityVector,unityVector,'k-')

% for igor

% tuning curves
identifier = ['SNmean1c',num2str(A)] ; 
ForIgor.(identifier) = SN1_mean ; 

identifier = ['SNmean2c',num2str(A)] ; 
ForIgor.(identifier) = SN2_mean ; 

% identifier = ['SNstd1c',num2str(A)] ; 
% ForIgor.(identifier) = SN1_std ; 
% 
% identifier = ['SNstd2c',num2str(A)] ; 
% ForIgor.(identifier) = SN2_std ; 

identifier = ['SNvar1c',num2str(A)] ; 
ForIgor.(identifier) = SN1_var ; 

identifier = ['SNvar2c',num2str(A)] ; 
ForIgor.(identifier) = SN2_var ; 


% % measured noise correlation coeficients
% for a=1:numUniqueBars ; % for each bar
%     identifier = ['snResLinCoef',num2str(BarAngles_CA_unique(a)),id,num2str(A)] ; 
%     ForIgor.(identifier) = SN_res_lincoef(a) ; 
% 
%     identifier = ['stResLinCoef',num2str(BarAngles_CA_unique(a)),id,num2str(A)] ;
%     ForIgor.(identifier) = ST_res_lincoef(a) ;
% end
%         
% identifier = ['snResLinCoefall',id,num2str(A)] ;
% ForIgor.(identifier) = SN_res_lincoef ;
% 
% identifier = ['stResLinCoefall',id,num2str(A)] ;
% ForIgor.(identifier) = ST_res_lincoef ;
% 
% % correlation histograms 
% identifier = ['snResLinCoefHist',id,num2str(A)] ;
% ForIgor.(identifier) = hist(SN_res_lincoef,[-1:.1:1]) ;
% 
% identifier = ['stResLinCoefHist',id,num2str(A)] ;
% ForIgor.(identifier) = hist(ST_res_lincoef,[-1:.1:1]) ;
% 
% % st differences
% identifier = ['stDiffMeanall',id,num2str(A)] ;
% ForIgor.(identifier) =  abs(STdiff_mean) ;
% 
% % bar angles
% identifier = ['BarAngles',num2str(A)] ;
% ForIgor.(identifier) =  BarAngles_CA_unique' ;
% 
% identifier = ['vectorAngleDiff',id,num2str(A)] ;
% ForIgor.(identifier) = vectorAngleDiff ;



end



