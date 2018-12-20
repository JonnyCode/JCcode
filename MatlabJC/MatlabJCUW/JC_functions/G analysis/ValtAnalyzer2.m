function ForIgor = ValtAnalyzer2(Input,Parameters,id,id2, A) ; 
% this function will get cc +/- Inh block and correct the cc for
% contamination
% id should be inh block trials and id2 should be control trials to be
% corrected

% get data
try
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);
catch
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);
end

for section=1:2 ; % for each section (id or id2)

    if section==2 ; %
       id = id2 ;       
    end

    epochs = str2num(Input(A).(id)) ;
    round = 0 ;
    for a = 1:3:length(epochs) ; % for each spike epoch
        round = round +1 ;

        [dataExc(round,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data

        [dataInh(round,:), error] = ITCReadEpoch(epochs(a+1), 0, fp) ;    % get data

        [dataAltV(round,:), error] = ITCReadEpoch(epochs(a+2), 0, fp) ;    % get data

        [lightCommand(round,:), error] = ITCReadEpochStm(epochs(a+2), 0,fp); % get light
        [voltageCommand(round,:), error] = ITCReadEpochStm(epochs(a+2), 1,fp); % get voltage command

        [SI(round), error] = ITCGetSamplingInterval(epochs(a+2), fp); % get sampling interval
        SI(round) = SI(round) * 1e-6; % Sampling interval in sec
    end

    if Input(A).ITC18flag == 1 ;
        SI = SI*1.25 ;
    end

    lightCommand = lightCommand(:,1:length(dataExc)) ;
    lightCommand(lightCommand<0)=0 ;

    voltageCommand = voltageCommand(:,1:length(dataExc)) ;

    cyclepnts = 100 ; % number of sample points in a cycle, pnts between leaving hold1 and returning (Also gets rid of first cycle)    
    FirstAltPnt = (cyclepnts/2)+1 ; % first sample point you want to plot after begining of step from alternation 
    LastAltPnt = FirstAltPnt ;  % last "                                                                    "

    [prePnts, error] = ITCGetStmPrePts(epochs(1), 0, 0, fp) ; % points collected beyond which data is not worth analyzing
    [postPnts, error] = ITCGetStmTailPts(epochs(1), 0, 0, fp) ;

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

    % low pass filter to remove electrical crap
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
        GExc_hpf(a,:) = dataExc_hpf(a,:)/-61 - mean(dataExc_hpf(a,prePnts:end-postPnts)/-61) ; % get conductance from stable currrents
        GInh_hpf(a,:) = dataInh_hpf(a,:)/61 - mean(dataInh_hpf(a,prePnts:end-postPnts)/61) ; 
        GAlt_ExcInt_hpf(a,:) = Alt_ExcInt_hpf(a,:)/-61 - mean(Alt_ExcInt_hpf(a,prePnts:end-postPnts)/-61) ; % get conductance from alt current
        GAlt_InhInt_hpf(a,:) = Alt_InhInt_hpf(a,:)/61 - mean(Alt_InhInt_hpf(a,prePnts:end-postPnts)/61) ; %#ok<*AGROW>

    end

    offsetExc = min(min([GExc_hpf(:,prePnts:end-postPnts),GAlt_ExcInt_hpf(:,prePnts:end-postPnts)])) ;
    offsetInh = min(min([GInh_hpf(:,prePnts:end-postPnts),GAlt_InhInt_hpf(:,prePnts:end-postPnts)])) ;
    % ofset G (assumes all g have same mean and 1 min) 
    for a = 1:size(dataExc,1) ; % for each trial
        GExc_hpf(a,:) = GExc_hpf(a,:) - offsetExc ; % offsets
        GInh_hpf(a,:) = GInh_hpf(a,:) - offsetInh ; 
        GAlt_ExcInt_hpf(a,:) = GAlt_ExcInt_hpf(a,:) - offsetExc ; 
        GAlt_InhInt_hpf(a,:) = GAlt_InhInt_hpf(a,:) - offsetInh ;

    end

    % get mean data
    GAlt_ExcInt_hpf_Mean = mean(GAlt_ExcInt_hpf) ;
    GAlt_InhInt_hpf_Mean = mean(GAlt_InhInt_hpf) ;

    GExc_hpf_Mean = mean(GExc_hpf) ;
    GInh_hpf_Mean = mean(GInh_hpf) ;

    % get residuals
    for a= 1:size(GExc_hpf,1) ;
        GExc_hpf_Res(a,:) =  GExc_hpf(a,:) - GExc_hpf_Mean ;
        GInh_hpf_Res(a,:) =  GInh_hpf(a,:) - GInh_hpf_Mean ;
        GAlt_ExcInt_hpf_Res(a,:) = GAlt_ExcInt_hpf(a,:) - GAlt_ExcInt_hpf_Mean ;
        GAlt_InhInt_hpf_Res(a,:) = GAlt_InhInt_hpf(a,:) - GAlt_InhInt_hpf_Mean ;
    end

    % get cross correlations (avoid first 3sec because of oscillations from highpass filtering) 
    start = floor(3/SI(1)) ; % 3seconds in points
    for a= 1:size(GExc_hpf,1) ;
        ccAltPre(a,:) = xcov(GAlt_ExcInt_hpf(a,start:prePnts),GAlt_InhInt_hpf(a,start:prePnts)) ; %cross corr of alternating g exc and inh prepoints 
        ccPre(a,:) = xcov(GExc_hpf(a,start:prePnts),GInh_hpf(a,start:prePnts)) ; % cross corr of single hold g exc and inh prepoints

        ccAltRes(a,:) = xcov(GAlt_ExcInt_hpf_Res(a,prePnts:end-postPnts),GAlt_InhInt_hpf_Res(a,prePnts:end-postPnts)) ; % cross corr of alternating g exc and inh stimulus residuals
        ccRes(a,:) = xcov(GExc_hpf_Res(a,prePnts:end-postPnts),GInh_hpf_Res(a,prePnts:end-postPnts)) ; % cross corr of single hold g exc and inh stimulus residuals    

        acExcPre(a,:) = xcov(GExc_hpf(a,start:prePnts)) ; % auttocorrelation of prestim exc single hold
        acInhPre(a,:) = xcov(GInh_hpf(a,start:prePnts)) ; % auttocorrelation of prestim inh single hold

        acAltExcRes(a,:) = xcov(GAlt_ExcInt_hpf_Res(a,prePnts:end-postPnts)) ; % auttocorrelation of residual exc alt hold
        acAltInhRes(a,:) = xcov(GAlt_InhInt_hpf_Res(a,prePnts:end-postPnts)) ; % auttocorrelation of residual inh alt hold

        ccAltRes_Coef(a,:) = xcov(GAlt_ExcInt_hpf_Res(a,prePnts:end-postPnts),GAlt_InhInt_hpf_Res(a,prePnts:end-postPnts),'coef') ;
    end

    ccAlt_Mean = xcov(GAlt_ExcInt_hpf_Mean(1,prePnts:end-postPnts),GAlt_InhInt_hpf_Mean(1,prePnts:end-postPnts)) ; % cross corr of alternating g exc and inh stimulus mean
    cc_Mean = xcov(GExc_hpf_Mean(1,prePnts:end-postPnts),GInh_hpf_Mean(1,prePnts:end-postPnts)) ; % cross corr of single hold g exc and inh stimulus mean

    % mean cross corr 
    ccAltPre_Mean = mean(ccAltPre) ;
    ccPre_Mean = mean(ccPre) ;

    ccAltRes_Mean = mean(ccAltRes) ;
    ccRes_Mean = mean(ccRes) ;

    acExcPre_Mean = mean(acExcPre) ;
    acInhPre_Mean = mean(acInhPre) ;

    acAltExcRes_Mean = mean(acAltExcRes) ;
    acAltInhRes_Mean = mean(acAltInhRes) ;

    acAltExcRes_Mean_peak = max(acAltExcRes_Mean) ;
    acAltInhRes_Mean_peak = max(acAltInhRes_Mean) ;

    ccAltRes_Coef_Mean = mean(ccAltRes_Coef) ; % as calculated throughout paper
    ccAltRes_Coef_SEM = std(ccAltRes_Coef)/sqrt(size(ccAltRes_Coef,1)) ;
    
    % save stuff from inh block analysis
    if section==1 ; % if its inh block trial
        acAltExcRes_Mean_peak_INH = acAltExcRes_Mean_peak ;
        acAltInhRes_Mean_peak_INH = acAltInhRes_Mean_peak ;
        
        ccAltRes_INH = ccAltRes ;
        ccAltRes_Mean_INH = ccAltRes_Mean ;
        
        clearvars -except acAltExcRes_Mean_peak_INH acAltInhRes_Mean_peak_INH ccAltRes_INH ccAltRes_Mean_INH Input A id id2 fp section
    
    else
        d = sqrt(acAltInhRes_Mean_peak_INH/acAltExcRes_Mean_peak) ;
        
        CorrelationDenominator = sqrt(((acAltInhRes_Mean_peak - (d^2).*acAltExcRes_Mean_peak + 2*d.*(CCpeakFinder(ccAltRes_Mean) - CCpeakFinder(ccAltRes_Mean_INH)))./(1-d).^2).*acAltExcRes_Mean_peak) ; % denominator corrected
        
        ccAltRes_Mean_Normalized = mean(ccAltRes/sqrt(acAltExcRes_Mean_peak*acAltInhRes_Mean_peak)) ; % normalized as normal (like coef) but denominator is not per trial (as it normaly is)
        
        ccAltRes_Mean_NormalizedCorrected = (ccAltRes_Mean/(1-d))/CorrelationDenominator ; % normalized by total variance of unconamitated g
        
        ccAltRes_INH_NormalizedCorrected = (ccAltRes_INH/(1-d))/CorrelationDenominator ; % normalized by total variance of unconamitated g
        ccAltRes_Mean_INH_NormalizedCorrected = mean(ccAltRes_INH_NormalizedCorrected) ; 
        ccAltRes_SEM_INH_NormalizedCorrected = std(ccAltRes_INH_NormalizedCorrected)/sqrt(size(ccAltRes_INH_NormalizedCorrected,1)) ;
        
        ccAltRes_Mean_TrueCoef = ccAltRes_Mean_NormalizedCorrected - ccAltRes_Mean_INH_NormalizedCorrected ; % noise correlation without contamination 

        ccAltRes_Normalized = ccAltRes/sqrt(acAltExcRes_Mean_peak*acAltInhRes_Mean_peak) ; % normalized as normal (like coef)
        ccAltRes_SEM_Normalized = std(ccAltRes_Normalized)/sqrt(size(ccAltRes_Normalized,1)) ; % sem
        
        fractionUnderestimate = (max(ccAltRes_Mean_TrueCoef)-max(ccAltRes_Coef_Mean))/max(ccAltRes_Mean_TrueCoef) ;
        
    end
    
end

%time vectors for cc
time_cc = [SI(1)*([1:length(ccAltRes)] - (length(ccAltRes)+1)/2)] ;

figure
subplot(1,2,1)
hold on
plot(time_cc,ccAltRes_Mean_TrueCoef,'k--')
plot(time_cc,ccAltRes_Coef_Mean,'k')
plot(time_cc,ccAltRes_Mean_Normalized,'c--')

subplot(1,2,2)
plot(time_cc,ccAltRes_Mean_TrueCoef,'k--')
hold on
plot(time_cc,ccAltRes_Mean_NormalizedCorrected,'k')
plot(time_cc,ccAltRes_Mean_INH_NormalizedCorrected,'r')

% %forIgor
% identifier = ['timeCC',num2str(A)] ;
% ForIgor.(identifier) = time_cc ;
% 
% 
% identifier = ['ccAltResMean',num2str(A)] ; % noise correlations as in paper
% ForIgor.(identifier) = ccAltRes_Coef_Mean ;
% 
% identifier = ['ccAltResSEM',num2str(A)] ;
% ForIgor.(identifier) = ccAltRes_Coef_SEM ;
% 
identifier = ['ccAltResPeak',num2str(A)] ;
ForIgor.(identifier) = CCpeakFinder(ccAltRes_Coef_Mean) ;
% 
% 
% 
% identifier = ['ccAltResInhMean',num2str(A)] ; % noise correlations contamination
% ForIgor.(identifier) = ccAltRes_Mean_INH_NormalizedCorrected ;
% 
% identifier = ['ccAltResInhSEM',num2str(A)] ;
% ForIgor.(identifier) = ccAltRes_SEM_INH_NormalizedCorrected ;
% 
identifier = ['ccAltResInhPeak',num2str(A)] ;
ForIgor.(identifier) = CCpeakFinder(ccAltRes_Mean_INH_NormalizedCorrected) ;
% 
% 
% 
% identifier = ['ccAltResCorrectedMean',num2str(A)] ; % noise correlations true after correction
% ForIgor.(identifier) = ccAltRes_Mean_TrueCoef ;
% 
identifier = ['ccAltResCorrectedPeak',num2str(A)] ;
ForIgor.(identifier) = CCpeakFinder(ccAltRes_Mean_TrueCoef) ;



identifier = ['fractionUnderestimate',num2str(A)] ;
ForIgor.(identifier) = fractionUnderestimate ;





