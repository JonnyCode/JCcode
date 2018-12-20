function ForIgor = pairAltDS2(Input,Parameters,id,A) ; 

% this function is modified from pairAltDS.m and will analyze data from pairs of ds cells recroded in alt
% voltage exp.  I have added estimate based on tuning curves and simple
% stim dependent gain model and removed converging corr estimates based on
% pair recordings.

% JC 6/24/11 
% JC 7/15/11 EDITED to use stim points for length of stim not data length - pre and post pnts as done previously
% JC 9/6/11 EDITED to correct for skipped triggers and other stimulus issues.

residualOption = Parameters.residualOption ;

% get data
[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);


epochs = str2num(Input(A).(id)) ;
round = 0 ;
for a = 1:3:length(epochs) ; % for each spike epoch
    round = round +1 ;

    [voltageCommand1(round,:), error] = ITCReadEpochStm(epochs(a+2), 3,fp) ; % get voltage command
    [voltageCommand2(round,:), error] = ITCReadEpochStm(epochs(a+2), 4,fp) ; 
    
    [dataAltV1(round,:), error] = ITCReadEpoch(epochs(a+2), 0, fp) ;    % get data
    [dataAltV2(round,:), error] = ITCReadEpoch(epochs(a+2), 1, fp) ;
    
    
    if voltageCommand1(round,1)==min(voltageCommand1(round,1:500)) ;
        [dataExc1(round,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    %#ok<*AGROW> % get data
        [dataInh1(round,:), error] = ITCReadEpoch(epochs(a+1), 0, fp) ;    % get data
    else
        [dataExc1(round,:), error] = ITCReadEpoch(epochs(a+1), 0, fp) ;    %#ok<*AGROW> % get data
        [dataInh1(round,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data
    end
        
    if voltageCommand2(round,1)==min(voltageCommand2(round,1:500)) ;
        [dataExc2(round,:), error] = ITCReadEpoch(epochs(a), 1, fp) ;    %#ok<*AGROW> % get data
        [dataInh2(round,:), error] = ITCReadEpoch(epochs(a+1), 1, fp) ;    % get data
    else
        [dataExc2(round,:), error] = ITCReadEpoch(epochs(a+1), 1, fp) ;    %#ok<*AGROW> % get data
        [dataInh2(round,:), error] = ITCReadEpoch(epochs(a), 1, fp) ;    % get data
    end      

    
    [SI(round), error] = ITCGetSamplingInterval(epochs(a+2), fp); % get sampling interval
    SI(round) = SI(round) * 1e-6; % Sampling interval in sec
end

voltageCommand1 = voltageCommand1(:,1:length(dataExc1)) ;
voltageCommand2 = voltageCommand2(:,1:length(dataExc1)) ;

if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end

SpatialStimParams = hdf5load(['~/Data/mouse/',Input(A).cellname,'_spatial.h5']) ; % load spatial stim params

for a = 1:length(epochs) ;
    [Monitor(a,:), error] = ITCReadEpoch(epochs(a), 2 , fp) ;    % get monitor data
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
LastAltPnt = FirstAltPnt ;  % last "                                                                    "

samplerate = 1/SI(1) ; % Hz at which data was collected
time = [SI(1):SI(1):SI(1)*length(dataExc1(1,:))] ; % time vector in seconds

% WINDOWING ALTERNATING V DATA
Alt_Exc1 = NaN(size(dataAltV1)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold1
Alt_Inh1 = NaN(size(dataAltV1)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold2

Alt_Exc2 = NaN(size(dataAltV2)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold1
Alt_Inh2 = NaN(size(dataAltV2)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold2

for a = FirstAltPnt:LastAltPnt ;        % for each data point per cycle we want to plot  
    if voltageCommand1(1,1)==min(voltageCommand1(1,1:500)) ; 
        Alt_Exc1(:,[a:cyclepnts:end]) = dataAltV1(:,[a:cyclepnts:end]) ; 
        Alt_Inh1(:,[a+cyclepnts/2:cyclepnts:end]) = dataAltV1(:,[a+cyclepnts/2:cyclepnts:end]) ; 
    else
        Alt_Inh1(:,[a:cyclepnts:end]) = dataAltV1(:,[a:cyclepnts:end]) ; 
        Alt_Exc1(:,[a+cyclepnts/2:cyclepnts:end]) = dataAltV1(:,[a+cyclepnts/2:cyclepnts:end]) ;  
    end

    if voltageCommand2(1,1)==min(voltageCommand2(1,1:500)) ; 
        Alt_Exc2(:,[a:cyclepnts:end]) = dataAltV2(:,[a:cyclepnts:end]) ; 
        Alt_Inh2(:,[a+cyclepnts/2:cyclepnts:end]) = dataAltV2(:,[a+cyclepnts/2:cyclepnts:end]) ; 
    else
        Alt_Inh2(:,[a:cyclepnts:end]) = dataAltV2(:,[a:cyclepnts:end]) ; 
        Alt_Exc2(:,[a+cyclepnts/2:cyclepnts:end]) = dataAltV2(:,[a+cyclepnts/2:cyclepnts:end]) ;
    end
end 

numTrials = length(epochs)/3 ;
% interpolate alternating current data
for a = 1:numTrials ; % for each trial
    b = find(isnan(Alt_Exc1(a,:)) == 0) ;     % find all the indices that are not nans
    Alt_ExcInt1(a,:) = interp1(b,Alt_Exc1(a,b),[1:length(Alt_Exc1)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b
    
    b = find(isnan(Alt_Exc2(a,:)) == 0) ;     % find all the indices that are not nans
    Alt_ExcInt2(a,:) = interp1(b,Alt_Exc2(a,b),[1:length(Alt_Exc2)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b
    
    b = find(isnan(Alt_Inh1(a,:)) == 0) ;     % find all the indices that are not nans
    Alt_InhInt1(a,:) = interp1(b,Alt_Inh1(a,b),[1:length(Alt_Inh1)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b
    
    b = find(isnan(Alt_Inh2(a,:)) == 0) ;     % find all the indices that are not nans
    Alt_InhInt2(a,:) = interp1(b,Alt_Inh2(a,b),[1:length(Alt_Inh2)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b
end

% low pass filter and remove electrical crap
dataExc_lpf1 = lowPassFilter(dataExc1, samplerate, 5000) ; %(signal,samplerate,cutoff frequ (hz))
dataInh_lpf1 = lowPassFilter(dataInh1, samplerate, 5000) ; 
Alt_ExcInt_lpf1 = lowPassFilter(Alt_ExcInt1, samplerate, 5000) ;
Alt_InhInt_lpf1 = lowPassFilter(Alt_InhInt1, samplerate, 5000) ;

dataExc_lpf2 = lowPassFilter(dataExc2, samplerate, 5000) ; %(signal,samplerate,cutoff frequ (hz))
dataInh_lpf2 = lowPassFilter(dataInh2, samplerate, 5000) ; 
Alt_ExcInt_lpf2 = lowPassFilter(Alt_ExcInt2, samplerate, 5000) ;
Alt_InhInt_lpf2 = lowPassFilter(Alt_InhInt2, samplerate, 5000) ;


% high pass filter to remove slow drift
% greg S. sent this to me to help implement a butterworth and avoid ringing
F=1 ; % filter cuttoff
Wn = F*SI(1); %normalized frequency cutoff
[z, p, k] = butter(1,Wn,'high'); %
[sos,g]=zp2sos(z,p,k); 
myfilt=dfilt.df2sos(sos,g);

dataExc_hpf1 = filter(myfilt,dataExc_lpf1')'; % filter implementation
dataInh_hpf1 = filter(myfilt,dataInh_lpf1')'; 
Alt_ExcInt_hpf1 = filter(myfilt,Alt_ExcInt_lpf1')'; 
Alt_InhInt_hpf1 = filter(myfilt,Alt_InhInt_lpf1')';

dataExc_hpf2 = filter(myfilt,dataExc_lpf2')'; % filter implementation
dataInh_hpf2 = filter(myfilt,dataInh_lpf2')'; 
Alt_ExcInt_hpf2 = filter(myfilt,Alt_ExcInt_lpf2')'; 
Alt_InhInt_hpf2 = filter(myfilt,Alt_InhInt_lpf2')';

% change to conductances and subtract off means 
for a = 1:numTrials ; % for each trial
    GExc_hpf1(a,:) = dataExc_hpf1(a,:)/-61 - mean(dataExc_hpf1(a,prePnts(a):prePnts(a)+stimPnts(a))/-61) ; % get conductance from stable currrents
    GInh_hpf1(a,:) = dataInh_hpf1(a,:)/61 - mean(dataInh_hpf1(a,prePnts(a):prePnts(a)+stimPnts(a))/61) ; 
    GAlt_ExcInt_hpf1(a,:) = Alt_ExcInt_hpf1(a,:)/-61 - mean(Alt_ExcInt_hpf1(a,prePnts(a):prePnts(a)+stimPnts(a))/-61) ; % get conductance from alt current
    GAlt_InhInt_hpf1(a,:) = Alt_InhInt_hpf1(a,:)/61 - mean(Alt_InhInt_hpf1(a,prePnts(a):prePnts(a)+stimPnts(a))/61) ; %#ok<*AGROW>
    
    GExc_hpf2(a,:) = dataExc_hpf2(a,:)/-61 - mean(dataExc_hpf2(a,prePnts(a):prePnts(a)+stimPnts(a))/-61) ; % get conductance from stable currrents
    GInh_hpf2(a,:) = dataInh_hpf2(a,:)/61 - mean(dataInh_hpf2(a,prePnts(a):prePnts(a)+stimPnts(a))/61) ; 
    GAlt_ExcInt_hpf2(a,:) = Alt_ExcInt_hpf2(a,:)/-61 - mean(Alt_ExcInt_hpf2(a,prePnts(a):prePnts(a)+stimPnts(a))/-61) ; % get conductance from alt current
    GAlt_InhInt_hpf2(a,:) = Alt_InhInt_hpf2(a,:)/61 - mean(Alt_InhInt_hpf2(a,prePnts(a):prePnts(a)+stimPnts(a))/61) ; %#ok<*AGROW>
end

offsetAltExc1 = min(min(GAlt_ExcInt_hpf1(:,prePnts(a):prePnts(a)+stimPnts(a)))) ;
offsetExc1 = min(min(GExc_hpf1(:,prePnts(a):prePnts(a)+stimPnts(a)))) ;

offsetAltInh1 = min(min(GAlt_InhInt_hpf1(:,prePnts(a):prePnts(a)+stimPnts(a)))) ;
offsetInh1 = min(min(GInh_hpf1(:,prePnts(a):prePnts(a)+stimPnts(a)))) ;

offsetAltExc2 = min(min(GAlt_ExcInt_hpf2(:,prePnts(a):prePnts(a)+stimPnts(a)))) ;
offsetExc2 = min(min(GExc_hpf2(:,prePnts(a):prePnts(a)+stimPnts(a)))) ;

offsetAltInh2 = min(min(GAlt_InhInt_hpf2(:,prePnts(a):prePnts(a)+stimPnts(a)))) ;
offsetInh2 = min(min(GInh_hpf2(:,prePnts(a):prePnts(a)+stimPnts(a)))) ;

% ofset G (assumes all g have same mean and 1 min) 
for a = 1:numTrials ; % for each trial
    GExc_hpf1(a,:) = GExc_hpf1(a,:) - offsetExc1 ; % offsets
    GInh_hpf1(a,:) = GInh_hpf1(a,:) - offsetInh1 ; 
    GAlt_ExcInt_hpf1(a,:) = GAlt_ExcInt_hpf1(a,:) - offsetAltExc1 ; 
    GAlt_InhInt_hpf1(a,:) = GAlt_InhInt_hpf1(a,:) - offsetAltInh1 ;
    
    GExc_hpf2(a,:) = GExc_hpf2(a,:) - offsetExc2 ; % offsets
    GInh_hpf2(a,:) = GInh_hpf2(a,:) - offsetInh2 ; 
    GAlt_ExcInt_hpf2(a,:) = GAlt_ExcInt_hpf2(a,:) - offsetAltExc2 ; 
    GAlt_InhInt_hpf2(a,:) = GAlt_InhInt_hpf2(a,:) - offsetAltInh2 ;
end

% figure % example alternating I exctraction
% subplot(2,1,1)
% plot(time,dataAltV1(3,:))
% hold on
% plot(time,Alt_Exc1(3,:),'go')
% plot(time,Alt_Inh1(3,:),'ro')
% plot(time,Alt_ExcInt1(3,:),'g-')
% plot(time,Alt_InhInt1(3,:),'r-')
% xlabel('time (seconds)')
% ylabel('current (pA)')
% title('alternating extraction')
% 
% subplot(2,1,2)
% plot(time,dataAltV2(3,:))
% hold on
% plot(time,Alt_Exc2(3,:),'go')
% plot(time,Alt_Inh2(3,:),'ro')
% plot(time,Alt_ExcInt2(3,:),'g-')
% plot(time,Alt_InhInt2(3,:),'r-')
% xlabel('time (seconds)')
% ylabel('current (pA)')
% title('alternating extraction')

% to save memmory clear unused variables
clearvars -except GExc_hpf1 GInh_hpf1 GAlt_ExcInt_hpf1 GAlt_InhInt_hpf1...
    GExc_hpf2 GInh_hpf2 GAlt_ExcInt_hpf2 GAlt_InhInt_hpf2...
    prePnts postPnts BarPnts OnPnts BarAngles NumBars SI Monitor A Input fp...
    frameRate SpatialStimParams id time numTrials epochs residualOption stimPnts

UniqueBars = unique(BarAngles(:)) ; % unique bar angles
NumUniqueBars = length(UniqueBars) ; % number of unique bar angles

for a=1:NumUniqueBars ;
    gExc1{a} = [] ;
    gInh1{a} = [] ;
    gAlt_Exc1{a} = [] ;
    gAlt_Inh1{a} = [] ;

    gExc2{a} = [] ;
    gInh2{a} = [] ; 
    gAlt_Exc2{a} = [] ;
    gAlt_Inh2{a} = [] ;
end

% cut up and arrange data in array by bar direction 
r=0;
[Ba,i] = sort(BarAngles,2) ; % sort each row of bar angle in accending order
for b= 1:3:size(BarAngles,1) ; % on each trial
    r=r+1 ;
    for a = 1:NumBars(b) ; % for each bar
        Ubar = find(UniqueBars==Ba(b,a)) ;
        
        gExc1{Ubar} = [gExc1{Ubar}; GExc_hpf1(r,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*i(b,a))] ;
        gInh1{Ubar} = [gInh1{Ubar}; GInh_hpf1(r,prePnts(b+1)+BarPnts(b+1)*(i(b+1,a)-1):prePnts(b+1)+BarPnts(b+1)*i(b+1,a))] ;
        gAlt_Exc1{Ubar} = [gAlt_Exc1{Ubar}; GAlt_ExcInt_hpf1(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a))] ;
        gAlt_Inh1{Ubar} = [gAlt_Inh1{Ubar}; GAlt_InhInt_hpf1(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a))] ;
    
        gExc2{Ubar} = [gExc2{Ubar}; GExc_hpf2(r,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*i(b,a))] ;
        gInh2{Ubar} = [gInh2{Ubar}; GInh_hpf2(r,prePnts(b+1)+BarPnts(b+1)*(i(b+1,a)-1):prePnts(b+1)+BarPnts(b+1)*i(b+1,a))] ; 
        gAlt_Exc2{Ubar} = [gAlt_Exc2{Ubar}; GAlt_ExcInt_hpf2(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a))] ;
        gAlt_Inh2{Ubar} = [gAlt_Inh2{Ubar}; GAlt_InhInt_hpf2(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a))] ;    
    end
end

numBarRepeats = size(gExc1{1},1) ;

% get mean data
for a = 1:NumUniqueBars ;
    gExc_Mean1{a} = mean(gExc1{a}) ;
    gInh_Mean1{a} = mean(gInh1{a}) ;
    gAlt_Exc_Mean1{a} = mean(gAlt_Exc1{a}) ;
    gAlt_Inh_Mean1{a} = mean(gAlt_Inh1{a}) ;
    
    gExc_Mean2{a} = mean(gExc2{a}) ;
    gInh_Mean2{a} = mean(gInh2{a}) ;
    gAlt_Exc_Mean2{a} = mean(gAlt_Exc2{a}) ;
    gAlt_Inh_Mean2{a} = mean(gAlt_Inh2{a}) ;
end

% get residuals
for a = 1:NumUniqueBars ;
    if residualOption == 0 ;
        gExc_Res1{a} =  gExc1{a} - repmat(gExc_Mean1{a},numBarRepeats,1) ;
        gInh_Res1{a} =  gInh1{a} - repmat(gInh_Mean1{a},numBarRepeats,1) ;
        gAlt_Exc_Res1{a} = gAlt_Exc1{a} - repmat(gAlt_Exc_Mean1{a},numBarRepeats,1) ;
        gAlt_Inh_Res1{a} = gAlt_Inh1{a} - repmat(gAlt_Inh_Mean1{a},numBarRepeats,1) ;

        gExc_Res2{a} =  gExc2{a} - repmat(gExc_Mean2{a},numBarRepeats,1) ;
        gInh_Res2{a} =  gInh2{a} - repmat(gInh_Mean2{a},numBarRepeats,1) ;
        gAlt_Exc_Res2{a} = gAlt_Exc2{a} - repmat(gAlt_Exc_Mean2{a},numBarRepeats,1) ;
        gAlt_Inh_Res2{a} = gAlt_Inh2{a} - repmat(gAlt_Inh_Mean2{a},numBarRepeats,1) ;
    
    elseif residualOption == 1 ;
        for b = 2:numBarRepeats-1 ;
            gExc_Res1{a}(b-1,:) =  gExc1{a}(b,:) - (gExc1{a}(b-1,:)+gExc1{a}(b+1,:))/2 ;
            gInh_Res1{a}(b-1,:) =  gInh1{a}(b,:) - (gInh1{a}(b-1,:)+gInh1{a}(b+1,:))/2 ;
            gAlt_Exc_Res1{a}(b-1,:) =  gAlt_Exc1{a}(b,:) - (gAlt_Exc1{a}(b-1,:)+gAlt_Exc1{a}(b+1,:))/2 ;
            gAlt_Inh_Res1{a}(b-1,:) =  gAlt_Inh1{a}(b,:) - (gAlt_Inh1{a}(b-1,:)+gAlt_Inh1{a}(b+1,:))/2 ;

            gExc_Res2{a}(b-1,:) =  gExc2{a}(b,:) - (gExc2{a}(b-1,:)+gExc2{a}(b+1,:))/2 ;
            gInh_Res2{a}(b-1,:) =  gInh2{a}(b,:) - (gInh2{a}(b-1,:)+gInh2{a}(b+1,:))/2 ;
            gAlt_Exc_Res2{a}(b-1,:) =  gAlt_Exc2{a}(b,:) - (gAlt_Exc2{a}(b-1,:)+gAlt_Exc2{a}(b+1,:))/2 ;
            gAlt_Inh_Res2{a}(b-1,:) =  gAlt_Inh2{a}(b,:) - (gAlt_Inh2{a}(b-1,:)+gAlt_Inh2{a}(b+1,:))/2 ;
        end
    end    
end
                
% g peaks for tunning curves 
for a= 1:NumUniqueBars ;
    gExc_Mean_max1(a) = max(mean(gExc1{a})) ;
    gInh_Mean_max1(a) = max(mean(gInh1{a})) ;
    gAlt_Exc_Mean_max1(a) = max(mean(gAlt_Exc1{a})) ;
    gAlt_Inh_Mean_max1(a) = max(mean(gAlt_Inh1{a})) ;
    
    gExc_Mean_max2(a) = max(mean(gExc2{a})) ;
    gInh_Mean_max2(a) = max(mean(gInh2{a})) ;
    gAlt_Exc_Mean_max2(a) = max(mean(gAlt_Exc2{a})) ;
    gAlt_Inh_Mean_max2(a) = max(mean(gAlt_Inh2{a})) ;
end
    
%time vectors for cc
time_cc = [SI(1)*([1:2*length(gAlt_Exc1{a})-1] - length(gAlt_Exc1{a}))] ;
   
% get cross correlations 
for a= 1:NumUniqueBars ; 
    
    for b = 1:numBarRepeats-residualOption*2;    
        % coefs
        ccRes1{a}(b,:) = xcorr(gExc_Res1{a}(b,:),gInh_Res1{a}(b,:),'coef') ; % cc of res exc1/inh1
        ccAltRes1{a}(b,:) = xcorr(gAlt_Exc_Res1{a}(b,:),gAlt_Inh_Res1{a}(b,:),'coef') ;

        ccRes2{a}(b,:) = xcorr(gExc_Res2{a}(b,:),gInh_Res2{a}(b,:),'coef') ; % cc of res exc2/inh2
        ccAltRes2{a}(b,:) = xcorr(gAlt_Exc_Res2{a}(b,:),gAlt_Inh_Res2{a}(b,:),'coef') ;
 
        ccExcRes{a}(b,:) = xcorr(gExc_Res1{a}(b,:),gExc_Res2{a}(b,:),'coef') ; % cc of res exc1/exc2
        ccAltExcRes{a}(b,:) = xcorr(gAlt_Exc_Res1{a}(b,:),gAlt_Exc_Res2{a}(b,:),'coef') ; 
        
        ccInhRes{a}(b,:) = xcorr(gInh_Res1{a}(b,:),gInh_Res2{a}(b,:),'coef') ; % cc of res inh1/inh2 
        ccAltInhRes{a}(b,:) = xcorr(gAlt_Inh_Res1{a}(b,:),gAlt_Inh_Res2{a}(b,:),'coef') ; 
        
        ccExcInhRes{a}(b,:) = xcorr(gExc_Res1{a}(b,:),gInh_Res2{a}(b,:),'coef') ; % cc of res exc1/inh2 
        ccAltExcInhRes{a}(b,:) = xcorr(gAlt_Exc_Res1{a}(b,:),gAlt_Inh_Res2{a}(b,:),'coef') ;        
        
        ccInhExcRes{a}(b,:) = xcorr(gExc_Res2{a}(b,:),gInh_Res1{a}(b,:),'coef') ; % cc of res exc2/inh1 
        ccAltInhExcRes{a}(b,:) = xcorr(gAlt_Exc_Res2{a}(b,:),gAlt_Inh_Res1{a}(b,:),'coef') ;                  
       
        % corrs
        cRes1{a}(b,:) = xcorr(gExc_Res1{a}(b,:),gInh_Res1{a}(b,:),'biased') ; % cc of res exc1/inh1
        cAltRes1{a}(b,:) = xcorr(gAlt_Exc_Res1{a}(b,:),gAlt_Inh_Res1{a}(b,:),'biased') ;

        cRes2{a}(b,:) = xcorr(gExc_Res2{a}(b,:),gInh_Res2{a}(b,:),'biased') ; % cc of res exc2/inh2
        cAltRes2{a}(b,:) = xcorr(gAlt_Exc_Res2{a}(b,:),gAlt_Inh_Res2{a}(b,:),'biased') ;
 
        cExcRes{a}(b,:) = xcorr(gExc_Res1{a}(b,:),gExc_Res2{a}(b,:),'biased') ; % cc of res exc1/exc2
        cAltExcRes{a}(b,:) = xcorr(gAlt_Exc_Res1{a}(b,:),gAlt_Exc_Res2{a}(b,:),'biased') ; 
        
        cInhRes{a}(b,:) = xcorr(gInh_Res1{a}(b,:),gInh_Res2{a}(b,:),'biased') ; % cc of res inh1/inh2 
        cAltInhRes{a}(b,:) = xcorr(gAlt_Inh_Res1{a}(b,:),gAlt_Inh_Res2{a}(b,:),'biased') ; 
        
        cExcInhRes{a}(b,:) = xcorr(gExc_Res1{a}(b,:),gInh_Res2{a}(b,:),'biased') ; % cc of res exc1/inh2 
        cAltExcInhRes{a}(b,:) = xcorr(gAlt_Exc_Res1{a}(b,:),gAlt_Inh_Res2{a}(b,:),'biased') ;        
        
        cInhExcRes{a}(b,:) = xcorr(gExc_Res2{a}(b,:),gInh_Res1{a}(b,:),'biased') ; % cc of res exc2/inh1 
        cAltInhExcRes{a}(b,:) = xcorr(gAlt_Exc_Res2{a}(b,:),gAlt_Inh_Res1{a}(b,:),'biased') ;      
        
    end
    
    % mean cross corr coefs
    ccResMean1(a,:) = mean(ccRes1{a}) ;
    ccAltResMean1(a,:) = mean(ccAltRes1{a}) ;

    ccResMean2(a,:) = mean(ccRes2{a}) ;
    ccAltResMean2(a,:) = mean(ccAltRes2{a}) ;    
    
    ccExcResMean(a,:) = mean(ccExcRes{a}) ;
    ccAltExcResMean(a,:) = mean(ccAltExcRes{a}) ;
    
    ccInhResMean(a,:) = mean(ccInhRes{a}) ;
    ccAltInhResMean(a,:) = mean(ccAltInhRes{a}) ;
    
    ccExcInhResMean(a,:) = mean(ccExcInhRes{a}) ;
    ccAltExcInhResMean(a,:) = mean(ccAltExcInhRes{a}) ;        

    ccInhExcResMean(a,:) = mean(ccInhExcRes{a}) ; 
    ccAltInhExcResMean(a,:) = mean(ccAltInhExcRes{a}) ; 

    % mean cross corr coefs
    cAltResMean1(a,:) = mean(cAltRes1{a}) ;
    cAltResMean2(a,:) = mean(cAltRes2{a}) ;    
    cAltExcResMean(a,:) = mean(cAltExcRes{a}) ;
    cAltInhResMean(a,:) = mean(cAltInhRes{a}) ;
    cAltExcInhResMean(a,:) = mean(cAltExcInhRes{a}) ;        
    cAltInhExcResMean(a,:) = mean(cAltInhExcRes{a}) ; 
    
    % cc peaks
    ccAltResPeak1(a) = CCpeakFinder(ccAltResMean1(a,:)) ;
    ccAltResPeak2(a) = CCpeakFinder(ccAltResMean2(a,:)) ;
    ccAltExcResMean_peak(a) = CCpeakFinder(ccAltExcResMean(a,:)) ;
    ccAltInhResMean_peak(a) = CCpeakFinder(ccAltInhResMean(a,:)) ;
    ccAltExcInhResMean_peak(a) = CCpeakFinder(ccAltExcInhResMean(a,:)) ;        
    ccAltInhExcResMean_peak(a) = CCpeakFinder(ccAltInhExcResMean(a,:)) ; 
    
    % c peaks
    cAltResMean1_peak(a) = CCpeakFinder(cAltResMean1(a,:)) ;
    cAltResMean2_peak(a) = CCpeakFinder(cAltResMean2(a,:)) ;    
    cAltExcResMean_peak(a) = CCpeakFinder(cAltExcResMean(a,:)) ;
    cAltInhResMean_peak(a) = CCpeakFinder(cAltInhResMean(a,:)) ;
    cAltExcInhResMean_peak(a) = CCpeakFinder(cAltExcInhResMean(a,:)) ;        
    cAltInhExcResMean_peak(a) = CCpeakFinder(cAltInhExcResMean(a,:)) ;  
    
    % estimating variance of diveging noise assuming simple model
    gainProduct_exc1exc2(a) = sqrt(max(xcov(gAlt_Exc_Mean1{a},'biased')))*sqrt(max(xcov(gAlt_Exc_Mean2{a},'biased'))) ;
    gainProduct_inh1inh2(a) = sqrt(max(xcov(gAlt_Inh_Mean1{a},'biased')))*sqrt(max(xcov(gAlt_Inh_Mean2{a},'biased'))) ; 
    gainProduct_exc1inh2(a) = sqrt(max(xcov(gAlt_Exc_Mean1{a},'biased')))*sqrt(max(xcov(gAlt_Inh_Mean2{a},'biased'))) ;  
    gainProduct_exc2inh1(a) = sqrt(max(xcov(gAlt_Exc_Mean2{a},'biased')))*sqrt(max(xcov(gAlt_Inh_Mean1{a},'biased'))) ; 
    
end

% noise corr coef peaks aligned by peak of tunigng curve
peakPlace = 4 ; % indicy of peak to which you should align

ccTuningX = [0:7]*45 ;
ccTuningX = circshift(ccTuningX,[0,peakPlace-1]) ;

peakPlaceShift = peakPlace - find(ccAltExcResMean_peak == max(ccAltExcResMean_peak)) ;
peakExcVector = circshift(ccAltExcResMean_peak,[0,peakPlaceShift]) ;

peakInhVectorExcAlgn = circshift(ccAltInhResMean_peak,[0,peakPlaceShift]) ;
peakExcInhVectorExcAlgn = circshift(ccAltExcInhResMean_peak,[0,peakPlaceShift]) ;
peakInhExcVectorExcAlgn = circshift(ccAltInhExcResMean_peak,[0,peakPlaceShift]) ;

peakPlaceShift = peakPlace - find(ccAltInhResMean_peak == max(ccAltInhResMean_peak)) ;
peakInhVector = circshift(ccAltInhResMean_peak,[0,peakPlaceShift]) ;

peakPlaceShift = peakPlace - find(ccAltExcInhResMean_peak == max(ccAltExcInhResMean_peak)) ;
peakExcInhVector = circshift(ccAltExcInhResMean_peak,[0,peakPlaceShift]) ;

peakPlaceShift = peakPlace - find(ccAltInhExcResMean_peak == max(ccAltInhExcResMean_peak)) ;
peakInhExcVector = circshift(ccAltInhExcResMean_peak,[0,peakPlaceShift]) ;

% gain products linear fits coefs
[temp,p] = corrcoef(gainProduct_exc1exc2,cAltExcResMean_peak) ;
gainProduct_exc1exc2_EstimateQuality_linCoef = temp(1,2) ; 
gainProduct_exc1exc2_EstimateQuality_linCoef_pValue = p(1,2) ; 

[temp,p] = corrcoef(gainProduct_inh1inh2,cAltInhResMean_peak) ;
gainProduct_inh1inh2_EstimateQuality_linCoef = temp(1,2) ; 
gainProduct_inh1inh2_EstimateQuality_linCoef_pValue = p(1,2) ; 

[temp,p] = corrcoef(gainProduct_exc1inh2,cAltExcInhResMean_peak) ;
gainProduct_exc1inh2_EstimateQuality_linCoef = temp(1,2) ; 
gainProduct_exc1inh2_EstimateQuality_linCoef_pValue = p(1,2) ; 

[temp,p] = corrcoef(gainProduct_exc2inh1,cAltInhExcResMean_peak) ;
gainProduct_exc2inh1_EstimateQuality_linCoef = temp(1,2) ; 
gainProduct_exc2inh1_EstimateQuality_linCoef_pValue = p(1,2) ; 

% gain product range
gainProduct_exc1exc2_range = range(gainProduct_exc1exc2) ;  
gainProduct_inh1inh2_range = range(gainProduct_inh1inh2) ; 
gainProduct_exc1inh2_range = range(gainProduct_exc1inh2) ; 
gainProduct_exc2inh1_range = range(gainProduct_exc2inh1) ; 

% linear fit parameters

fitCoefs = polyfit(gainProduct_exc1exc2,cAltExcResMean_peak,1) ;
gainProduct_exc1exc2_Fit_slope = fitCoefs(1) ; 
gainProduct_exc1exc2_Fit_yint = fitCoefs(2) ;
gainProduct_exc1exc2_line = fitCoefs(1)*gainProduct_exc1exc2 + fitCoefs(2) ;

fitCoefs = polyfit(gainProduct_inh1inh2,cAltInhResMean_peak,1) ;
gainProduct_inh1inh2_Fit_slope = fitCoefs(1) ;
gainProduct_inh1inh2_Fit_yint = fitCoefs(2) ; 
gainProduct_inh1inh2_line = fitCoefs(1)*gainProduct_inh1inh2 + fitCoefs(2) ;

fitCoefs = polyfit(gainProduct_exc1inh2,cAltExcInhResMean_peak,1) ;
gainProduct_exc1inh2_Fit_slope = fitCoefs(1) ; 
gainProduct_exc1inh2_Fit_yint = fitCoefs(2) ; 
gainProduct_exc1inh2_line = fitCoefs(1)*gainProduct_exc1inh2 + fitCoefs(2) ;

fitCoefs = polyfit(gainProduct_exc2inh1,cAltInhExcResMean_peak,1) ;
gainProduct_exc2inh1_Fit_slope = fitCoefs(1) ; 
gainProduct_exc2inh1_Fit_yint = fitCoefs(2) ; 
gainProduct_exc2inh1_line = fitCoefs(1)*gainProduct_exc2inh1 + fitCoefs(2) ;

% error of fit coefs (according to Taylor - "An intro. to error analysis" book)
gainProduct_exc1exc2_Fit_SquaredError = (cAltExcResMean_peak - gainProduct_exc1exc2_line).^2 ; 
gainProduct_exc1exc2_Fit_EstimatedStd = sqrt((1/(NumUniqueBars-2))*sum(gainProduct_exc1exc2_Fit_SquaredError)) ;
gainProduct_exc1exc2_Fit_slope_std = gainProduct_exc1exc2_Fit_EstimatedStd *sqrt(NumUniqueBars/(NumUniqueBars*sum(gainProduct_exc1exc2.^2) - sum(gainProduct_exc1exc2)^2)) ;
gainProduct_exc1exc2_Fit_yint_std = gainProduct_exc1exc2_Fit_EstimatedStd *sqrt((sum(gainProduct_exc1exc2.^2))/(NumUniqueBars*sum(gainProduct_exc1exc2.^2) - sum(gainProduct_exc1exc2)^2)) ;

gainProduct_inh1inh2_Fit_SquaredError = (cAltExcResMean_peak - gainProduct_inh1inh2_line).^2 ; 
gainProduct_inh1inh2_Fit_EstimatedStd = sqrt((1/(NumUniqueBars-2))*sum(gainProduct_inh1inh2_Fit_SquaredError)) ;
gainProduct_inh1inh2_Fit_slope_std = gainProduct_inh1inh2_Fit_EstimatedStd *sqrt(NumUniqueBars/(NumUniqueBars*sum(gainProduct_inh1inh2.^2) - sum(gainProduct_inh1inh2)^2)) ;
gainProduct_inh1inh2_Fit_yint_std = gainProduct_inh1inh2_Fit_EstimatedStd *sqrt((sum(gainProduct_inh1inh2.^2))/(NumUniqueBars*sum(gainProduct_inh1inh2.^2) - sum(gainProduct_inh1inh2)^2)) ;

gainProduct_exc1inh2_Fit_SquaredError = (cAltExcResMean_peak - gainProduct_exc1inh2_line).^2 ; 
gainProduct_exc1inh2_Fit_EstimatedStd = sqrt((1/(NumUniqueBars-2))*sum(gainProduct_exc1inh2_Fit_SquaredError)) ;
gainProduct_exc1inh2_Fit_slope_std = gainProduct_exc1inh2_Fit_EstimatedStd *sqrt(NumUniqueBars/(NumUniqueBars*sum(gainProduct_exc1inh2.^2) - sum(gainProduct_exc1inh2)^2)) ;
gainProduct_exc1inh2_Fit_yint_std = gainProduct_exc1inh2_Fit_EstimatedStd *sqrt((sum(gainProduct_exc1inh2.^2))/(NumUniqueBars*sum(gainProduct_exc1inh2.^2) - sum(gainProduct_exc1inh2)^2)) ;

gainProduct_exc2inh1_Fit_SquaredError = (cAltExcResMean_peak - gainProduct_exc2inh1_line).^2 ; 
gainProduct_exc2inh1_Fit_EstimatedStd = sqrt((1/(NumUniqueBars-2))*sum(gainProduct_exc2inh1_Fit_SquaredError)) ;
gainProduct_exc2inh1_Fit_slope_std = gainProduct_exc2inh1_Fit_EstimatedStd *sqrt(NumUniqueBars/(NumUniqueBars*sum(gainProduct_exc2inh1.^2) - sum(gainProduct_exc2inh1)^2)) ;
gainProduct_exc2inh1_Fit_yint_std = gainProduct_exc2inh1_Fit_EstimatedStd *sqrt((sum(gainProduct_exc2inh1.^2))/(NumUniqueBars*sum(gainProduct_exc2inh1.^2) - sum(gainProduct_exc2inh1)^2)) ;


% concatinated fit parameters
FitParamtersSlope = [gainProduct_exc1exc2_Fit_slope, gainProduct_inh1inh2_Fit_slope, gainProduct_exc1inh2_Fit_slope, gainProduct_exc2inh1_Fit_slope] ;
FitParamtersSlope_std = [gainProduct_exc1exc2_Fit_slope_std, gainProduct_inh1inh2_Fit_slope_std, gainProduct_exc1inh2_Fit_slope_std, gainProduct_exc2inh1_Fit_slope_std] ;

FitParamtersYint = [gainProduct_exc1exc2_Fit_yint, gainProduct_inh1inh2_Fit_yint, gainProduct_exc1inh2_Fit_yint, gainProduct_exc2inh1_Fit_yint] ;
FitParamtersYint_std = [gainProduct_exc1exc2_Fit_yint_std, gainProduct_inh1inh2_Fit_yint_std, gainProduct_exc1inh2_Fit_yint_std, gainProduct_exc2inh1_Fit_yint_std] ;

% concatinating data
ExcG1 = cell2mat(gAlt_Exc1) ;
InhG1 = cell2mat(gAlt_Inh1) ;
ExcG2 = cell2mat(gAlt_Exc2) ;
InhG2 = cell2mat(gAlt_Inh2) ;


% prep G for dynamic clamp

maxExc1 = max(max(GAlt_ExcInt_hpf1(:,prePnts(1):prePnts(1)+stimPnts(1)))) ;
maxExc2 = max(max(GAlt_ExcInt_hpf2(:,prePnts(1):prePnts(1)+stimPnts(1)))) ;

maxInh1 = max(max(GAlt_InhInt_hpf1(:,prePnts(1):prePnts(1)+stimPnts(1)))) ;
maxInh2 = max(max(GAlt_InhInt_hpf2(:,prePnts(1):prePnts(1)+stimPnts(1)))) ;

Gtime = [1:length(ExcG1)]*SI(1) ;


% % figures
% figure % monitor check
% plot(time,Monitor)
% 
% 
% figure % conductance check
% subplot(4,1,1)
% plot(Gtime,ExcG1)
% 
% subplot(4,1,2)
% plot(Gtime,InhG1)
% 
% subplot(4,1,3)
% plot(Gtime,ExcG2)
% 
% subplot(4,1,4)
% plot(Gtime,InhG2)


% % 
% % 
% figure % conductance check
% for a=1:numTrials ; % on each trial
%     subplot(4,1,1)
%     plot(time,GExc_hpf1(a,:))
%     hold on
%     plot(time,GInh_hpf1(a,:),'r')
%     hold off
%     
%     subplot(4,1,2)
%     plot(time,GExc_hpf2(a,:))
%     hold on
%     plot(time,GInh_hpf2(a,:),'r')
%     hold off
%     
%     subplot(4,1,3)
%     plot(time,GAlt_ExcInt_hpf1(a,:))
%     hold on
%     plot(time,GAlt_InhInt_hpf1(a,:),'r')
%     hold off
%     
%     subplot(4,1,4)
%     plot(time,GAlt_ExcInt_hpf2(a,:))
%     hold on
%     plot(time,GAlt_InhInt_hpf2(a,:),'r')
%     hold off
%     
%     text(0,.9,num2str(epochs((a-1)*3+1)),'units','norm')
%     pause
% end
% 
% 
figure
for a=1:NumUniqueBars ;
    subplot(ceil(NumUniqueBars(1)/2),2,a)
    plot(time_cc,ccAltResMean1(a,:),'k')
    hold on
    plot(time_cc,ccResMean1(a,:),'g')
end

figure
for a=1:NumUniqueBars ;
    subplot(ceil(NumUniqueBars(1)/2),2,a)
    plot(time_cc,ccAltResMean2(a,:),'k')
    hold on
    plot(time_cc,ccResMean2(a,:),'g')
end

figure
for a=1:NumUniqueBars ;
    subplot(ceil(NumUniqueBars(1)/2),2,a)
    plot(time_cc,ccAltExcResMean(a,:),'b')
    hold on
    plot(time_cc,ccAltInhResMean(a,:),'r')
end

figure
for a=1:NumUniqueBars ;
    subplot(ceil(NumUniqueBars(1)/2),2,a)
    plot(time_cc,ccExcResMean(a,:),'b')
    hold on
    plot(time_cc,ccInhResMean(a,:),'r')
end

figure
for a=1:NumUniqueBars ;
    subplot(ceil(NumUniqueBars(1)/2),2,a)
    plot(time_cc,ccAltInhExcResMean(a,:),'m')
    hold on
    plot(time_cc,ccAltExcInhResMean(a,:),'c')
end
% 
% figure
% subplot(1,2,1)
% plot(peakExcVector)
% hold on 
% plot(peakInhVector,'r')
% plot(peakExcInhVector,'c')
% plot(peakInhExcVector,'m')
% 
% subplot(1,2,2)
% plot(peakExcVector)
% hold on 
% plot(peakInhVectorExcAlgn,'r')
% plot(peakExcInhVectorExcAlgn,'c')
% plot(peakInhExcVectorExcAlgn,'m')

% 
% figure
% plot(lowerBound1)
% hold on
% plot(upperBound1)
% plot(Estimate,'b--')
% plot(ccAltResPeak1,'k')
% plot(ccAltResPeak2,'k--')

figure % simple gain product prediction
subplot(2,4,1:4)
plot(gainProduct_exc1exc2,cAltExcResMean_peak,'*b')
hold on
plot(gainProduct_inh1inh2,cAltInhResMean_peak,'*r')
plot(gainProduct_exc1inh2,cAltExcInhResMean_peak,'*c')
plot(gainProduct_exc2inh1,cAltInhExcResMean_peak,'*m')

plot(gainProduct_exc1exc2,gainProduct_exc1exc2_line,'b')
plot(gainProduct_inh1inh2,gainProduct_inh1inh2_line,'r')
plot(gainProduct_exc1inh2,gainProduct_exc1inh2_line,'c')
plot(gainProduct_exc2inh1,gainProduct_exc2inh1_line,'m')

xlabel('gain (auttocorr peak of mean response)')
ylabel('covariance (biased cross corr peak of residuals)')

subplot(2,4,5)
plot(1,gainProduct_exc1exc2_EstimateQuality_linCoef,'b*')
hold on
plot(1,gainProduct_inh1inh2_EstimateQuality_linCoef,'r*')
plot(1,gainProduct_exc1inh2_EstimateQuality_linCoef,'c*')
plot(1,gainProduct_exc2inh1_EstimateQuality_linCoef,'m*')
ylabel('corr coef')

subplot(2,4,6)
plot(1,gainProduct_exc1exc2_EstimateQuality_linCoef_pValue,'b*')
hold on
plot(1,gainProduct_inh1inh2_EstimateQuality_linCoef_pValue,'r*')
plot(1,gainProduct_exc1inh2_EstimateQuality_linCoef_pValue,'c*')
plot(1,gainProduct_exc2inh1_EstimateQuality_linCoef_pValue,'m*')
ylabel('corr coef p Value')

subplot(2,4,7)
errorbar(1,gainProduct_exc1exc2_Fit_slope,gainProduct_exc1exc2_Fit_slope_std,'b*')
hold on
errorbar(1,gainProduct_inh1inh2_Fit_slope,gainProduct_inh1inh2_Fit_slope_std,'r*')
errorbar(1,gainProduct_exc1inh2_Fit_slope,gainProduct_exc1inh2_Fit_slope_std,'c*')
errorbar(1,gainProduct_exc2inh1_Fit_slope,gainProduct_exc2inh1_Fit_slope_std,'m*')
ylabel('slope')

subplot(2,4,8)
errorbar(1,gainProduct_exc1exc2_Fit_yint,gainProduct_exc1exc2_Fit_yint_std,'b*')
hold on
errorbar(1,gainProduct_inh1inh2_Fit_yint,gainProduct_inh1inh2_Fit_yint_std,'r*')
errorbar(1,gainProduct_exc1inh2_Fit_yint,gainProduct_exc1inh2_Fit_yint_std,'c*')
errorbar(1,gainProduct_exc2inh1_Fit_yint,gainProduct_exc2inh1_Fit_yint_std,'m*')
ylabel('y intercept')

% 
% figure
% plot(gainProduct_exc1exc2_range,gainProduct_exc1exc2_EstimateQuality_linCoef,'*') ;  
% hold on
% plot(gainProduct_inh1inh2_range,gainProduct_inh1inh2_EstimateQuality_linCoef,'r*') ;
% plot(gainProduct_exc1inh2_range,gainProduct_exc1inh2_EstimateQuality_linCoef,'c*') ; 
% plot(gainProduct_exc2inh1_range,gainProduct_exc2inh1_EstimateQuality_linCoef,'m*') ; 
% 

% For igor

identifier = ['GtimeCell',num2str(A)] ;
ForIgor.(identifier) = Gtime ;

for b=1:numBarRepeats ;

%    % normalized for dynamic clamp
%     identifier = ['ExcG1t',num2str(b),'c',num2str(A)] ;
%     ForIgor.(identifier) = ExcG1(b,:)/maxExc1 ;
%     
%     identifier = ['InhG1t',num2str(b),'c',num2str(A)] ;
%     ForIgor.(identifier) = InhG1(b,:)/maxInh1 ;
%        
%     identifier = ['ExcG2t',num2str(b),'c',num2str(A)] ;
%     ForIgor.(identifier) = ExcG2(b,:)/maxExc2 ;
%        
%     identifier = ['InhG2t',num2str(b),'c',num2str(A)] ;
%     ForIgor.(identifier) = InhG2(b,:)/maxInh2 ;


    % conductances not normalized
    identifier = ['ExcG1t',num2str(b),'cell',num2str(A)] ;
    ForIgor.(identifier) = ExcG1(b,:) ;
    
    identifier = ['InhG1t',num2str(b),'cell',num2str(A)] ;
    ForIgor.(identifier) = InhG1(b,:) ;
       
    identifier = ['ExcG2t',num2str(b),'cell',num2str(A)] ;
    ForIgor.(identifier) = ExcG2(b,:) ;
       
    identifier = ['InhG2t',num2str(b),'cell',num2str(A)] ;
    ForIgor.(identifier) = InhG2(b,:) ;
end


% cross correlations
identifier = ['timecc',num2str(A)] ;
ForIgor.(identifier) = time_cc ;

for a=1:NumUniqueBars ;
    identifier = ['ccAltResMean1b',num2str(a),'c',num2str(A)] ;
    ForIgor.(identifier) = ccAltResMean1(a,:) ;

    identifier = ['ccResMean1b',num2str(a),'c',num2str(A)] ;
    ForIgor.(identifier) = ccResMean1(a,:) ;

    identifier = ['ccAltResMean2b',num2str(a),'c',num2str(A)] ;
    ForIgor.(identifier) = ccAltResMean2(a,:) ;

    identifier = ['ccResMean2b',num2str(a),'c',num2str(A)] ;
    ForIgor.(identifier) = ccResMean2(a,:) ;

    identifier = ['ccAltExcResMeanb',num2str(a),'c',num2str(A)] ;
    ForIgor.(identifier) = ccAltExcResMean(a,:) ;

    identifier = ['ccAltInhResMeanb',num2str(a),'c',num2str(A)] ;
    ForIgor.(identifier) = ccAltInhResMean(a,:) ;

    identifier = ['ccExcResMeanb',num2str(a),'c',num2str(A)] ;
    ForIgor.(identifier) = ccExcResMean(a,:) ;
    
    identifier = ['ccInhResMeanb',num2str(a),'c',num2str(A)] ;
    ForIgor.(identifier) = ccInhResMean(a,:) ;

    identifier = ['ccAltExcInhResMeanb',num2str(a),'c',num2str(A)] ;
    ForIgor.(identifier) = ccAltExcInhResMean(a,:) ;

    identifier = ['ccAltInhExcResMeanb',num2str(a),'c',num2str(A)] ;
    ForIgor.(identifier) = ccAltInhExcResMean(a,:) ;

    identifier = ['ccExcInhResMeanb',num2str(a),'c',num2str(A)] ;
    ForIgor.(identifier) = ccExcInhResMean(a,:) ;

    identifier = ['ccInhExcResMeanb',num2str(a),'c',num2str(A)] ;
    ForIgor.(identifier) = ccInhExcResMean(a,:) ;

end

% cross correlation peak vectors

identifier = ['ccTuningXCell',num2str(A)] ;
ForIgor.(identifier) = ccTuningX ;

identifier = ['ccTuningExcCell',num2str(A)] ;
ForIgor.(identifier) = peakExcVector ;

identifier = ['ccTuningInhCell',num2str(A)] ;
ForIgor.(identifier) = peakInhVector ;

identifier = ['ccTuningExcInhCell',num2str(A)] ;
ForIgor.(identifier) = peakExcInhVector ;

identifier = ['ccTuningInhExcCell',num2str(A)] ;
ForIgor.(identifier) = peakInhExcVector ;

identifier = ['ccTuningInhAlgnExcCell',num2str(A)] ;
ForIgor.(identifier) = peakInhVectorExcAlgn ;

identifier = ['ccTuningExcInhAlgnExcCell',num2str(A)] ;
ForIgor.(identifier) = peakExcInhVectorExcAlgn ;

identifier = ['ccTuningInhExcAlgnExcCell',num2str(A)] ;
ForIgor.(identifier) = peakInhExcVectorExcAlgn ;


% gain factor plots

identifier = ['gainPexcCell',num2str(A)] ;
ForIgor.(identifier) = gainProduct_exc1exc2 ;

identifier = ['gainPinhCell',num2str(A)] ;
ForIgor.(identifier) = gainProduct_inh1inh2 ;

identifier = ['gainPexcInhCell',num2str(A)] ;
ForIgor.(identifier) = gainProduct_exc1inh2 ;

identifier = ['gainPinhExcCell',num2str(A)] ;
ForIgor.(identifier) = gainProduct_exc2inh1 ;



identifier = ['cPeakExc',num2str(A)] ;
ForIgor.(identifier) = cAltExcResMean_peak ;

identifier = ['cPeakInh',num2str(A)] ;
ForIgor.(identifier) = cAltInhResMean_peak ;

identifier = ['cPeakExcInh',num2str(A)] ;
ForIgor.(identifier) = cAltExcInhResMean_peak ;

identifier = ['cPeakInhExc',num2str(A)] ;
ForIgor.(identifier) = cAltInhExcResMean_peak ;

% fit correlations and gainproduct range
identifier = ['gpExcLineCoefCell',num2str(A)] ;
ForIgor.(identifier) = gainProduct_exc1exc2_EstimateQuality_linCoef ;

identifier = ['gpInhLineCoefCell',num2str(A)] ;
ForIgor.(identifier) = gainProduct_inh1inh2_EstimateQuality_linCoef ;

identifier = ['gpExcInhLineCoefCell',num2str(A)] ;
ForIgor.(identifier) = gainProduct_exc1inh2_EstimateQuality_linCoef ;

identifier = ['gpInhExcLineCoefCell',num2str(A)] ;
ForIgor.(identifier) = gainProduct_exc2inh1_EstimateQuality_linCoef ;


identifier = ['gpExcRangeCell',num2str(A)] ;
ForIgor.(identifier) = gainProduct_exc1exc2_range ;

identifier = ['gpInhRangeCell',num2str(A)] ;
ForIgor.(identifier) = gainProduct_inh1inh2_range ;

identifier = ['gpExcInhRangeCell',num2str(A)] ;
ForIgor.(identifier) = gainProduct_exc1inh2_range ;

identifier = ['gpInhExcRangeCell',num2str(A)] ;
ForIgor.(identifier) = gainProduct_exc2inh1_range ;

% fit parameters

identifier = ['fitSlope',num2str(A)] ;
ForIgor.(identifier) = FitParamtersSlope ;

identifier = ['fitSlopeStd',num2str(A)] ;
ForIgor.(identifier) = FitParamtersSlope_std ;

identifier = ['fitYint',num2str(A)] ;
ForIgor.(identifier) = FitParamtersYint ;

identifier = ['fitYintStd',num2str(A)] ;
ForIgor.(identifier) = FitParamtersYint_std ;

