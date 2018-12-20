function ForIgor = ValtAnalyzerDS3(Input,Parameters,id,A)  
% this function will analyze alternate voltage exp from ds cells presented
% a moving bar

% JC 8/18/11 adapted from 'ValtAnlyzerDS.m' (not ValtAnalyzerDS2.m) to deal
% with data recorded from single cells but during two amp and one amp modes

% JC 9/6/11 removed all on and off specific code

if ~isempty(Input(A).channel) ; % if this is a two amp situation and a channel is specified
    dataSegment = Input(A).channel ; % Channel you want to analyze if two channel
    monitorSegment = 2 ;
else
    dataSegment = 0 ;
    monitorSegment = 1 ;
end

residualOption = Parameters.residualOption ;

% get data
[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
round = 0 ;
for a = 1:3:length(epochs) ; % for each spike epoch
    round = round +1 ;

    [dataExc(round,:), error] = ITCReadEpoch(epochs(a), dataSegment, fp) ;    %#ok<*AGROW> % get data

    [dataInh(round,:), error] = ITCReadEpoch(epochs(a+1), dataSegment, fp) ;    % get data
    
    [dataAltV(round,:), error] = ITCReadEpoch(epochs(a+2), dataSegment, fp) ;    % get data
    
    [voltageCommand(round,:), error] = ITCReadEpochStm(epochs(a+2), dataSegment,fp); % get voltage command
    
    [SI(round), error] = ITCGetSamplingInterval(epochs(a+2), fp); % get sampling interval
    SI(round) = SI(round) * 1e-6; % Sampling interval in sec
end

voltageCommand = voltageCommand(:,1:length(dataExc)) ;

if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end

SpatialStimParams = hdf5load(['~/Data/mouse/',Input(A).cellname,'_spatial.h5']) ; % load spatial stim params

for a = 1:length(epochs) ;
    [Monitor(a,:), error] = ITCReadEpoch(epochs(a), monitorSegment , fp) ;    % get monitor data
end

safeUpsweep = 2 ; % the upsweep at which you trust every epoch triggered properly (does not include entirely missed triggers)

FrameStartPnts = MonitorAnalysis2(Monitor, length(epochs), SI(1)) ;
for a = 1:length(epochs) ; 
    framerate_mean(a) = 60/mean(diff(FrameStartPnts{a}(safeUpsweep:end))) ; % frames/pnt
end
frameRate = mean(framerate_mean) ;

% correcting each epoch for triggering issues if neccessary
if ~isempty(Input(A).triggerMiss) ; % if any epochs missed a trigger 
    triggerMiss = Input(A).triggerMiss ; % number of upsweeps missed for each epoch
    if length(triggerMiss) == 1 ; % if only one number is given than all epochs missed the same number of triggers
        triggerMiss = ones(1,length(epochs))*triggerMiss ;
    end
else
    triggerMiss = zeros(1,length(epochs)) ;
end

numSkippedUpsweeps = triggerMiss + safeUpsweep - 1 ; % total number of frame monitor upsweeps that were skipped
numSkippedPnts = floor((numSkippedUpsweeps*60)/frameRate) ;

for a = 1:length(epochs) ;
    FirstSafePnt(a) = FrameStartPnts{a}(safeUpsweep) ;
end

% frames to pnts
for a = 1:length(epochs) ;
    StrucString = ['params_epoch',num2str(epochs(a))] ; 
    Struct = SpatialStimParams.(StrucString) ;

    BarAngles(a,:) = Struct.BarAngle ;

    prePnts(a) = floor(Struct.spatial_prepts/frameRate) ;
    stimPnts(a) = floor(Struct.spatial_stimpts/frameRate) ;
    postPnts(a) = floor(Struct.spatial_postpts/frameRate) ;
    
    % correct for triggering issues
    if numSkippedPnts(a)<=prePnts(a) ; % if triggered within the prepnts
        prePnts(a) = FirstSafePnt(a) + (prePnts(a) - numSkippedPnts(a)) ;
    else
        error('numSkippedPnts is more than prePnts')
    end
 
    NumBars(a) = length(BarAngles(a,:)) ; % number of bars shown per trial
    BarPnts(a) = floor(stimPnts(a)/NumBars(a)) ;
end

cyclepnts = 100 ; % number of sample points in a cycle, pnts between leaving hold1 and returning (Also gets rid of first cycle)    
FirstAltPnt = (cyclepnts/2)+1 ; % first sample point you want to plot after begining of step from alternation 
%FirstAltPnt = 30 ;
LastAltPnt = FirstAltPnt ;  % last "                                                                    "

samplerate = 1/SI(1) ; % Hz at which data was collected
time = [SI(1):SI(1):SI(1)*length(dataExc(1,:))] ; % time vector in seconds

% WINDOWING ALTERNATING V DATA
Alt_Exc = NaN(size(dataAltV)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold1
Alt_Inh = NaN(size(dataAltV)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold2_

for a = FirstAltPnt:LastAltPnt ;        % for each data point per cycle we want to plot
    Alt_Exc(:,[a:cyclepnts:end]) = dataAltV(:,[a:cyclepnts:end]) ; 
    Alt_Inh(:,[a+cyclepnts/2:cyclepnts:end]) = dataAltV(:,[a+cyclepnts/2:cyclepnts:end]) ; 
end 

% interpolate alternating current data
for a = 1:size(dataExc,1) ; % for each trial
    b = find(isnan(Alt_Exc(a,:)) == 0) ;     % find all the indices that are not nans
    Alt_ExcInt(a,:) = interp1(b,Alt_Exc(a,b),[1:length(Alt_Exc)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b
    
    b = find(isnan(Alt_Inh(a,:)) == 0) ;     % find all the indices that are not nans
    Alt_InhInt(a,:) = interp1(b,Alt_Inh(a,b),[1:length(Alt_Inh)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b
end

% low pass filter and remove electrical crap
dataExc_lpf = lowPassFilter(dataExc, samplerate, 5000) ; %(signal,samplerate,cutoff frequ (hz))
dataInh_lpf = lowPassFilter(dataInh, samplerate, 5000) ; 
Alt_ExcInt_lpf = lowPassFilter(Alt_ExcInt, samplerate, 5000) ;
Alt_InhInt_lpf = lowPassFilter(Alt_InhInt, samplerate, 5000) ;

% high pass filter to remove slow drift
% greg S. sent this to me to help implement a butterworth and avoid ringing
F=1 ; % filter cuttoff
Wn = F*SI(1); %normalized frequency cutoff
[z, p, k] = butter(1,Wn,'high'); %
[sos,g]=zp2sos(z,p,k); 
myfilt=dfilt.df2sos(sos,g);

dataExc_hpf = filter(myfilt,dataExc_lpf')'; % filter implementation
dataInh_hpf = filter(myfilt,dataInh_lpf')'; 
Alt_ExcInt_hpf = filter(myfilt,Alt_ExcInt_lpf')'; 
Alt_InhInt_hpf = filter(myfilt,Alt_InhInt_lpf')';

% change to conductances and subtract off means 
for a = 1:size(dataExc,1) ; % for each trial
    GExc_hpf(a,:) = dataExc_hpf(a,:)/-61 - mean(dataExc_hpf(a,prePnts:prePnts+stimPnts)/-61) ; % get conductance from stable currrents
    GInh_hpf(a,:) = dataInh_hpf(a,:)/61 - mean(dataInh_hpf(a,prePnts:prePnts+stimPnts)/61) ; 
    GAlt_ExcInt_hpf(a,:) = Alt_ExcInt_hpf(a,:)/-61 - mean(Alt_ExcInt_hpf(a,prePnts:prePnts+stimPnts)/-61) ; % get conductance from alt current
    GAlt_InhInt_hpf(a,:) = Alt_InhInt_hpf(a,:)/61 - mean(Alt_InhInt_hpf(a,prePnts:prePnts+stimPnts)/61) ; %#ok<*AGROW>

end

offsetExc = min(min([GExc_hpf(:,prePnts:prePnts+stimPnts),GAlt_ExcInt_hpf(:,prePnts:prePnts+stimPnts)])) ;
offsetInh = min(min([GInh_hpf(:,prePnts:prePnts+stimPnts),GAlt_InhInt_hpf(:,prePnts:prePnts+stimPnts)])) ;
% ofset G (assumes all g have same mean and 1 min) 
for a = 1:size(dataExc,1) ; % for each trial
    GExc_hpf(a,:) = GExc_hpf(a,:) - offsetExc ; % offsets
    GInh_hpf(a,:) = GInh_hpf(a,:) - offsetInh ; 
    GAlt_ExcInt_hpf(a,:) = GAlt_ExcInt_hpf(a,:) - offsetExc ; 
    GAlt_InhInt_hpf(a,:) = GAlt_InhInt_hpf(a,:) - offsetInh ;
end

% to save memmory clear unused variables
clearvars -except GExc_hpf GInh_hpf GAlt_ExcInt_hpf GAlt_InhInt_hpf...
    prePnts BarPnts BarAngles NumBars SI Monitor A Input fp...
    frameRate SpatialStimParams id id2 time residualOption data* voltageCommand Alt*

% cut up and arrange data in array by bar direction 
r=0;
[Ba,i] = sort(BarAngles,2) ; % sort each row of bar angle in accending order
for b= 1:3:size(BarAngles,1) ; % on each trial
    r=r+1 ;
    for  a = 1:NumBars(1) ; % for each bar
        
        gExc{a}(r,:) = GExc_hpf(r,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*i(b,a)) ;
        gInh{a}(r,:) = GInh_hpf(r,prePnts(b+1)+BarPnts(b+1)*(i(b+1,a)-1):prePnts(b+1)+BarPnts(b+1)*i(b+1,a)) ;
        gAlt_Exc{a}(r,:) = GAlt_ExcInt_hpf(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a)) ;
        gAlt_Inh{a}(r,:) = GAlt_InhInt_hpf(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a)) ;
    
    end
end

% get mean, meanamp, and mean time peak  data
for a = 1:NumBars(1) ;
    gExc_Mean{a} = mean(gExc{a}) ;
    gInh_Mean{a} = mean(gInh{a}) ;
    gAlt_Exc_Mean{a} = mean(gAlt_Exc{a}) ;
    gAlt_Inh_Mean{a} = mean(gAlt_Inh{a}) ;
end

% get residuals
if residualOption == 0 ; % added residual option 8/4/10
    % standard
    for a= 1:NumBars(1) ;
        gExc_Res{a} =  gExc{a} - repmat(gExc_Mean{a},size(gExc{a},1),1) ;
        gInh_Res{a} =  gInh{a} - repmat(gInh_Mean{a},size(gInh{a},1),1) ;
        gAlt_Exc_Res{a} = gAlt_Exc{a} - repmat(gAlt_Exc_Mean{a},size(gAlt_Exc{a},1),1) ;
        gAlt_Inh_Res{a} = gAlt_Inh{a} - repmat(gAlt_Inh_Mean{a},size(gAlt_Inh{a},1),1) ;
    end

elseif residualOption == 1 ;
    % nearest neighbor residuals
    for a= 1:NumBars(1) ;
        for b=2:size(gExc{a},1)-1 ; % for each possible residual trial
            gExc_Res{a}(b-1,:) =  gExc{a}(b,:) - (gExc{a}(b-1,:)+gExc{a}(b+1,:))/2 ;
            gInh_Res{a}(b-1,:) =  gInh{a}(b,:) - (gInh{a}(b-1,:)+gInh{a}(b+1,:))/2 ;
            gAlt_Exc_Res{a}(b-1,:) = gAlt_Exc{a}(b,:) - (gAlt_Exc{a}(b-1,:)+gAlt_Exc{a}(b+1,:))/2 ;
            gAlt_Inh_Res{a}(b-1,:) = gAlt_Inh{a}(b,:) - (gAlt_Inh{a}(b-1,:)+gAlt_Inh{a}(b+1,:))/2 ;

        end
        
    end
end
    
    
% g peaks for tunning curves 
for a= 1:NumBars(1) ;
    gExc_max_Mean(a) = mean(max(gExc{a},[],2)) ;
    gInh_max_Mean(a) = mean(max(gInh{a},[],2)) ;
    gAlt_Exc_max_Mean(a) = mean(max(gAlt_Exc{a},[],2)) ;
    gAlt_Inh_max_Mean(a) = mean(max(gAlt_Inh{a},[],2)) ;
        
    gExc_max_std(a) = std(max(gExc{a},[],2)) ;
    gInh_max_std(a) = std(max(gInh{a},[],2)) ;
    gAlt_Exc_max_std(a) = std(max(gAlt_Exc{a},[],2)) ;
    gAlt_Inh_max_std(a) = std(max(gAlt_Inh{a},[],2)) ;
    
end
    
%time vectors for cc
time_cc = [SI(1)*([1:2*length(gAlt_Exc{a})-1] - length(gAlt_Exc{a}))] ;

% get cross correlations 
for a= 1:NumBars(1) ;
    cc(a,:) = xcorr(gExc_Mean{a},gInh_Mean{a},'coef') ;
    ccAlt(a,:) = xcorr(gAlt_Exc_Mean{a},gAlt_Inh_Mean{a},'coef') ;
    
    gAlt_Inh_ResShuff{a} = circshift(gAlt_Inh_Res{a},[1,0]) ;

    for b = 1:size(gExc_Res{a},1) ;
    
        ccAltInd{a}(b,:) = xcorr(gAlt_Exc{a}(b,:),gAlt_Inh{a}(b,:),'coef') ; % cc of the ind alt (these may be different than the cc of the mean g above)
        ccRes{a}(b,:) = xcorr(gExc_Res{a}(b,:),gInh_Res{a}(b,:),'coef') ; % cc of res
        ccAltRes{a}(b,:) = xcorr(gAlt_Exc_Res{a}(b,:),gAlt_Inh_Res{a}(b,:),'coef') ;
        ccAltResShuff{a}(b,:) = xcorr(gAlt_Exc_Res{a}(b,:),gAlt_Inh_ResShuff{a}(b,:),'coef') ; % cc of res shuffled

    end
    % mean cross corr 
    ccAltIndMean(a,:) = mean(ccAltInd{a}) ;
    ccResMean(a,:) = mean(ccRes{a}) ;
    ccAltResMean(a,:) = mean(ccAltRes{a}) ;
    ccAltResShuffMean(a,:) = mean(ccAltResShuff{a}) ;
    
    % standard error of mean (sem) of cross correlations 
    ccRes_SEM(a,:) = std(ccRes{a})/sqrt(size(GExc_hpf,1)) ;
    ccAltRes_SEM(a,:) = std(ccAltRes{a})/sqrt(size(GExc_hpf,1)) ;
    ccAltResShuff_SEM(a,:) = std(ccAltResShuff{a})/sqrt(size(GExc_hpf,1)) ;    

    % cc peaks
    [ccRes_Peakabs,ccRes_Peaki] = max(abs(ccResMean(a,time_cc<.2 & time_cc>-.2))) ;
    ccRes_Peak(a) = ccResMean(a,ccRes_Peaki+find(time_cc>-.2,1)) ;

    [ccAltRes_Peakabs,ccAltRes_Peaki] = max(abs(ccAltResMean(a,time_cc<.2 & time_cc>-.2))) ;
    ccAltRes_Peak(a) = ccAltResMean(a,ccAltRes_Peaki+find(time_cc>-.2,1)) ;
    
%     % corr coef of meanamp
%     [cCoef,p] = corrcoef(gExc_ON_meanamp_Res{a},gInh_ON_meanamp_Res{a}) ;
%     corrCoef_meanamp(a) = cCoef(1,2) ;
%     p_meanamp(a) = p(1,2) ;
end
peakVector = [ccAltRes_Peak,nan(1,8-length(ccAltRes_Peak))] ; % make vector to align peaks
peakPlace = 4 ; % indicy of peak to which you should align
peakPlaceShift = peakPlace - find(ccAltRes_Peak == max(ccAltRes_Peak)) ;
peakVector = circshift(peakVector,[0,peakPlaceShift]) ;

% estimating covariance and variance of diveging noise assuming simple model
for a= 1:NumBars(1) ;
    % residual autocorrs
    for b = 1:size(gExc_Res{a},1) ;
        ccAltRes_forEstimate{a}(b,:) = xcorr(gAlt_Exc_Res{a}(b,:),gAlt_Inh_Res{a}(b,:),'biased') ;
        ccAltRes_excAc_forEstimate{a}(b,:) = xcorr(gAlt_Exc_Res{a}(b,:),'biased') ;
        ccAltRes_inhAc_forEstimate{a}(b,:) = xcorr(gAlt_Inh_Res{a}(b,:),'biased') ;
    end
    
    ccAltRes_forEstimateMean(a,:) =  mean(ccAltRes_forEstimate{a}) ;
    [temp_Peakabs,temp_Peaki] = max(abs(ccAltRes_forEstimateMean(a,time_cc<.2 & time_cc>-.2))) ;
    ccAltRes_forEstimate_Peak(a) = ccAltRes_forEstimateMean(a,temp_Peaki+find(time_cc>-.2,1)) ; % covariance peak
    
    acAltExcRes_forEstimate_Peak(a) = max(mean(ccAltRes_excAc_forEstimate{a})) ; % variance peaks
    acAltInhRes_forEstimate_Peak(a) = max(mean(ccAltRes_inhAc_forEstimate{a})) ;
    
    
    % conductance autocorrs
    ccAlt_excAc_forEstimate(a,:) = xcov(mean(gAlt_Exc{a}),'biased') ; % autto correlation of average response
    ccAlt_inhAc_forEstimate(a,:) = xcov(mean(gAlt_Inh{a}),'biased') ;

    acAltExc_forEstimate_Peak(a) = max(ccAlt_excAc_forEstimate(a,:)) ; % peaks of autocorr
    acAltInh_forEstimate_Peak(a) = max(ccAlt_inhAc_forEstimate(a,:)) ;

end

gainExc_squared = acAltExc_forEstimate_Peak ;
gainInh_squared = acAltInh_forEstimate_Peak ;
gainProduct = sqrt(gainExc_squared).*sqrt(gainInh_squared) ;

[temp,p] = corrcoef(gainProduct,ccAltRes_forEstimate_Peak) ;
covEstimateQuality_linCoef = temp(1,2) ; % ability to estimate cov from gain product
covEstimateQuality_pValue = p(1,2) ;

[temp,p] = corrcoef(gainExc_squared,acAltExcRes_forEstimate_Peak) ;
varExcEstimateQuality_linCoef = temp(1,2) ; % ability to estimate var from gain
varExcEstimateQuality_pValue = p(1,2) ;

[temp,p] = corrcoef(gainInh_squared,acAltInhRes_forEstimate_Peak) ;
varInhEstimateQuality_linCoef = temp(1,2) ; % ability to estimate var from gain
varInhEstimateQuality_pValue = p(1,2) ;

% range of gain factors
gainExc_squared_range = range(gainExc_squared) ;
gainInh_squared_range = range(gainInh_squared) ;
gainProduct_range = range(gainProduct) ;

% linear fits coefs

fitCoefs = polyfit(gainProduct,ccAltRes_forEstimate_Peak,1) ;
gainProductFit_slope = fitCoefs(1) ; % best fit common noise variance from covariance estimate
gainProductFit_yint = fitCoefs(2) ;
gainProductFit_line = gainProduct*gainProductFit_slope + gainProductFit_yint ;

fitCoefs = polyfit(gainExc_squared,acAltExcRes_forEstimate_Peak,1) ;
gainExc_squaredFit_slope = fitCoefs(1) ; % best fit common noise variance from variance estimate
gainExc_squaredFit_yint = fitCoefs(2) ; % best fit of independant exc noise variance from variance estimate
gainExc_squaredFit_line = gainExc_squared*gainExc_squaredFit_slope + gainExc_squaredFit_yint ;

fitCoefs = polyfit(gainInh_squared,acAltInhRes_forEstimate_Peak,1) ;
gainInh_squaredFit_slope = fitCoefs(1) ; % best fit common noise variance from variance estimate
gainInh_squaredFit_yint = fitCoefs(2) ; % best fit of independant exc noise variance from variance estimate
gainInh_squaredFit_line = gainInh_squared*gainInh_squaredFit_slope + gainInh_squaredFit_yint ;

% error of fit coefs (according to Taylor - "An intro. to error analysis" book)
gainProductFit_SquaredError = (ccAltRes_forEstimate_Peak - gainProductFit_line).^2 ; 
gainProductFit_EstimatedStd = sqrt((1/(NumBars(1)-2))*sum(gainProductFit_SquaredError)) ;
gainProductFit_slope_std = gainProductFit_EstimatedStd*sqrt(NumBars(1)/(NumBars(1)*sum(gainProduct.^2) - sum(gainProduct)^2)) ;
gainProductFit_yint_std = gainProductFit_EstimatedStd*sqrt((sum(gainProduct.^2))/(NumBars(1)*sum(gainProduct.^2) - sum(gainProduct)^2)) ;

gainExc_squaredFit_SquaredError = (ccAltRes_forEstimate_Peak - gainExc_squaredFit_line).^2 ; 
gainExc_squaredFit_EstimatedStd = sqrt((1/(NumBars(1)-2))*sum(gainExc_squaredFit_SquaredError)) ;
gainExc_squaredFit_slope_std = gainExc_squaredFit_EstimatedStd*sqrt(NumBars(1)/(NumBars(1)*sum(gainExc_squared.^2) - sum(gainExc_squared)^2)) ;
gainExc_squaredFit_yint_std = gainExc_squaredFit_EstimatedStd*sqrt((sum(gainExc_squared.^2))/(NumBars(1)*sum(gainExc_squared.^2) - sum(gainExc_squared)^2)) ;

gainInh_squaredFit_SquaredError = (ccAltRes_forEstimate_Peak - gainInh_squaredFit_line).^2 ; 
gainInh_squaredFit_EstimatedStd = sqrt((1/(NumBars(1)-2))*sum(gainInh_squaredFit_SquaredError)) ;
gainInh_squaredFit_slope_std = gainInh_squaredFit_EstimatedStd*sqrt(NumBars(1)/(NumBars(1)*sum(gainInh_squared.^2) - sum(gainInh_squared)^2)) ;
gainInh_squaredFit_yint_std = gainInh_squaredFit_EstimatedStd*sqrt((sum(gainInh_squared.^2))/(NumBars(1)*sum(gainInh_squared.^2) - sum(gainInh_squared)^2)) ;

FitParamtersSlope = [gainProductFit_slope, gainExc_squaredFit_slope, gainInh_squaredFit_slope] ;
FitParamtersSlope_std = [gainProductFit_slope_std, gainExc_squaredFit_slope_std, gainInh_squaredFit_slope_std] ;

FitParamtersYint = [gainProductFit_yint, gainExc_squaredFit_yint, gainInh_squaredFit_yint] ;
FitParamtersYint_std = [gainProductFit_yint_std, gainExc_squaredFit_yint_std, gainInh_squaredFit_yint_std] ;

% to save memmory clear unused variables
% clearvars -except Monitor gAlt_Exc_ON_Mean gExc_ON_Mean gAlt_Inh_ON_Mean...
%     gInh_ON_Mean gAlt_Exc_OFF_Mean gAlt_Inh_OFF_Mean gExc_OFF_Mean gInh_OFF_Mean...
%     time_cc ccAltResMean ccResMean ccAltResShuffMean ccRes_SEM ccAltRes_SEM ccAltResShuff_SEM...
%     time_ccON cc_ON ccAlt_ON time_ccOFF cc_OFF ccAlt_OFF ccAltRes_ON_Mean...
%     ccAltResShuff_ON_Mean ccAltResShuff_ON_Mean ccAltRes_OFF_Mean...
%     ccAltResShuff_OFF_Mean Ba gAlt_Exc_ON_max_Mean gAlt_Inh_ON_max_Mean...
%     gExc_ON_max_Mean gInh_ON_max_Mean gAlt_Exc_OFF_max_Mean gAlt_Inh_OFF_max_Mean...
%     gExc_OFF_max_Mean gInh_OFF_max_Mean NumBars SI A Input fp SpatialStimParams...
%     frameRate ccRes_ON_Mean ccRes_OFF_Mean gAlt_Exc gAlt_Inh gAlt_Exc_ON gExc_ON gAlt_Inh_ON...
%     gInh_ON gAlt_Exc_OFF gAlt_Inh_OFF gExc_OFF gInh_OFF gAlt_Exc_Res gAlt_Inh_Res...
%     id2 id ccAltInd_ON_Mean ccAltInd_OFF_Mean time ccRes ccAltRes_Peak ccRes_Peak...
%     gAlt_Exc_max_Mean gAlt_Inh_max_Mean gAlt_Exc_max_std gAlt_Inh_max_std...
%     residualOption gAlt_Exc_Mean gAlt_Inh_Mean gainProduct gainExc_squared gainInh_squared...
%     ccAltRes_forEstimate_Peak acAltExcRes_forEstimate_Peak acAltInhRes_forEstimate_Peak...
%     covEstimateQuality_linCoef varExcEstimateQuality_linCoef varInhEstimateQuality_linCoef...
%     CommonNoiseVar_estimate ccAltRes_Peak peakVector
%     

% idividual alternating conductances, residuals and res cross correlations
concat_gAlt_Exc = cell2mat(gAlt_Exc) ; % concatinate bar directions
concat_gAlt_Inh = cell2mat(gAlt_Inh) ;

concat_gAlt_ExcRes = cell2mat(gAlt_Exc_Res) ;
concat_gAlt_InhRes = cell2mat(gAlt_Inh_Res) ;

maxExc = max(max(concat_gAlt_Exc)) ;
maxInh = max(max(concat_gAlt_Inh)) ;

concat_gAlt_Exc_Mean = cell2mat(gAlt_Exc_Mean) ; 
concat_gAlt_Inh_Mean = cell2mat(gAlt_Inh_Mean) ;

concat_gExc_Mean = cell2mat(gExc_Mean) ;
concat_gInh_Mean = cell2mat(gInh_Mean) ;

% figures
% figure
% subplot(2,1,1)
% plot(concat_gExc_Mean)
% hold on
% plot(concat_gAlt_Exc_Mean,'r')
% 
% subplot(2,1,2)
% plot(concat_gInh_Mean)
% hold on
% plot(concat_gAlt_Inh_Mean,'r')
% title(num2str(A))

% figure
% subplot(3,1,1)
% plot([1:length(Monitor)],Monitor)
% xlabel('sample pnts')
% ylabel('Monitor reading')
% title('conductances')

% figure
% plot(time(1:length(concat_gAlt_Exc)),concat_gAlt_Exc) 
% 
% figure
% plot(time(1:length(concat_gAlt_Inh)),concat_gAlt_Inh) 
% 
figure  
for a=1:8
    subplot(1,8,a)
    plot(time_cc,ccAltResMean(a,:),'k')
    hold on
    plot(time_cc,ccResMean(a,:),'g')
    set(gca,'xlim',[-.2,.2],'ylim',[-.5,1])
end
% 
figure
subplot(2,4,1:4)
plot(gainProduct,ccAltRes_forEstimate_Peak,'k*')
hold on
plot(gainProduct,gainProductFit_line,'k-')

plot(gainExc_squared,acAltExcRes_forEstimate_Peak,'*')
plot(gainExc_squared,gainExc_squaredFit_line)

plot(gainInh_squared,acAltInhRes_forEstimate_Peak,'r*')
plot(gainInh_squared,gainInh_squaredFit_line,'r-')

subplot(2,4,5)
plot(1,covEstimateQuality_linCoef,'k*')
hold on
plot(1,varExcEstimateQuality_linCoef,'b*')
plot(1,varInhEstimateQuality_linCoef,'r*')
ylabel('corr coef')

subplot(2,4,6)
plot(1,covEstimateQuality_pValue,'k*')
hold on
plot(1,varExcEstimateQuality_pValue,'b*')
plot(1,varInhEstimateQuality_pValue,'r*')
ylabel('corr coef p Value')

subplot(2,4,7)
errorbar(1,gainProductFit_slope,gainProductFit_slope_std,'k*')
hold on
errorbar(1,gainExc_squaredFit_slope,gainExc_squaredFit_slope_std,'b*')
errorbar(1,gainInh_squaredFit_slope,gainInh_squaredFit_slope_std,'r*')
ylabel('slope')

subplot(2,4,8)
errorbar(1,gainProductFit_yint,gainProductFit_yint_std,'k*')
hold on
errorbar(1,gainExc_squaredFit_yint,gainExc_squaredFit_yint_std,'b*')
errorbar(1,gainInh_squaredFit_yint,gainInh_squaredFit_yint_std,'r*')
ylabel('y intercept')





% for igor 

% % alt voltage current example (for methods figure)
% identifier = ['altIexample1cell',num2str(A)] ;
% ForIgor.(identifier) = dataAltV(1,:) ;
% 
% identifier = ['altIexcExample1cell',num2str(A)] ;
% ForIgor.(identifier) = Alt_ExcInt(1,:) ;
% 
% identifier = ['altIinhExample1cell',num2str(A)] ;
% ForIgor.(identifier) = Alt_InhInt(1,:) ;
% 
% identifier = ['AltIVoltage','cell',num2str(A)] ;
% ForIgor.(identifier) = voltageCommand(1,:) ;
% 
% identifier = ['AltItime','cell',num2str(A)] ;
% ForIgor.(identifier) = time ;



%mean cross correlations for each bar direction
% for a=1:8 ;
%     identifier = ['ccAltResMean',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = ccAltResMean(a,:) ;
%     
%     identifier = ['ccResMean',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = ccResMean(a,:) ;
%     
%     
%     identifier = ['ccAltResSEM',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = ccAltRes_SEM(a,:) ;
%     
%     identifier = ['ccResSEM',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = ccRes_SEM(a,:) ;
%     
% end
% 
% identifier = ['time_cc','cell',num2str(A)] ;
% ForIgor.(identifier) = time_cc ;

% % mean conductances
% identifier = ['meanGaltExc','cell',num2str(A)] ;
% ForIgor.(identifier) = concat_gAlt_Exc_Mean ;
% 
% identifier = ['meanGaltInh','cell',num2str(A)] ;
% ForIgor.(identifier) = concat_gAlt_Inh_Mean ;
% 
% 
% identifier = ['meanGExc','cell',num2str(A)] ;
% ForIgor.(identifier) = concat_gExc_Mean ;
% 
% identifier = ['meanGInh','cell',num2str(A)] ;
% ForIgor.(identifier) = concat_gInh_Mean ;

% idividual alternating conductances, residuals and res cross correlations
% for a=1:size(concat_gAlt_Exc,1) ; 
%     
% %     % normalized for g clamp
% %     identifier = ['ExcG1t',num2str(a),'c',num2str(A)] ;
% %     ForIgor.(identifier) = concat_gAlt_Exc(a,:)/maxExc ;
% % 
% %     identifier = ['InhG1t',num2str(a),'c',num2str(A)] ;
% %     ForIgor.(identifier) = concat_gAlt_Inh(a,:)/maxInh ;
% 
%     
%     % unormalized
%     identifier = ['exGaltExc',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = concat_gAlt_Exc(a,:) ;
% 
%     identifier = ['exGaltInh',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = concat_gAlt_Inh(a,:) ;
% end

% for a=1:size(concat_gAlt_ExcRes,1) ; 
% 
%     identifier = ['exGaltExcRes',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = concat_gAlt_ExcRes(a,:) ;
%     
%     identifier = ['exGaltInhRes',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = concat_gAlt_InhRes(a,:) ;    
% end

% identifier = ['time','cell',num2str(A)] ;
% ForIgor.(identifier) = [1:length(concat_gAlt_Exc)]*SI(1) ;
 

% estimating noise correlation based on tuning curves and simple model
% identifier = ['covExcInh','cell',num2str(A)] ;
% ForIgor.(identifier) = ccAltRes_forEstimate_Peak ;
% 
% identifier = ['varExc','cell',num2str(A)] ;
% ForIgor.(identifier) = acAltExcRes_forEstimate_Peak ;
% 
% identifier = ['varInh','cell',num2str(A)] ;
% ForIgor.(identifier) = acAltInhRes_forEstimate_Peak ;
% 
% identifier = ['gainProduct','cell',num2str(A)] ;
% ForIgor.(identifier) = gainProduct ;
% 
% identifier = ['gainExcSquared','cell',num2str(A)] ;
% ForIgor.(identifier) = gainExc_squared ;
% 
% identifier = ['gainInhSquared','cell',num2str(A)] ;
% ForIgor.(identifier) = gainInh_squared ;


identifier = ['CovEstimateLinCorr','cell',num2str(A)] ;
ForIgor.(identifier) = covEstimateQuality_linCoef ;

identifier = ['varExcEstimateLinCorr','cell',num2str(A)] ;
ForIgor.(identifier) = varExcEstimateQuality_linCoef ;

identifier = ['varInhEstimateLinCorr','cell',num2str(A)] ;
ForIgor.(identifier) = varInhEstimateQuality_linCoef ;

identifier = ['CovEstimateLinCorrPV','cell',num2str(A)] ;
ForIgor.(identifier) = covEstimateQuality_pValue ;

identifier = ['varExcEstimateLinCorrPV','cell',num2str(A)] ;
ForIgor.(identifier) = varExcEstimateQuality_pValue ;

identifier = ['varInhEstimateLinCorrPV','cell',num2str(A)] ;
ForIgor.(identifier) = varInhEstimateQuality_pValue ;


identifier = ['fitParamsSlope','cell',num2str(A)] ;
ForIgor.(identifier) = FitParamtersSlope ;

identifier = ['fitParamsSlopeStd','cell',num2str(A)] ;
ForIgor.(identifier) = FitParamtersSlope_std ;

identifier = ['fitParamsYint','cell',num2str(A)] ;
ForIgor.(identifier) = FitParamtersYint ;

identifier = ['fitParamsYintStd','cell',num2str(A)] ;
ForIgor.(identifier) = FitParamtersYint_std ;


identifier = ['gainExcRange','cell',num2str(A)] ;
ForIgor.(identifier) = gainExc_squared_range ;

identifier = ['gainInhRange','cell',num2str(A)] ;
ForIgor.(identifier) = gainInh_squared_range ;

identifier = ['gainProductRange','cell',num2str(A)] ;
ForIgor.(identifier) = gainProduct_range ;
% 
% tuning curve of peak cc aligned by peak to a 1:8 vector peak at 4
identifier = ['ccPeakVector','cell',num2str(A)] ;
ForIgor.(identifier) = peakVector ;
% 

