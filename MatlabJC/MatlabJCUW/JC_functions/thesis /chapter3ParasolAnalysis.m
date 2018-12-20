function ForIgor = chapter3ParasolAnalysis(Input,A) ;

% thesis analysis for chapter 3 parasol data (cc of exc and inh)

% JC 8/22/11

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);
 
epochsExc = str2num(Input(A).Exc) ;
epochsInh = str2num(Input(A).Inh) ;
epochsCA = str2num(Input(A).CA) ;

for a = 1:length(epochsExc) ; % for each g epoch
    [dataExc(a,:),error] = ITCReadEpoch(epochsExc(a), 0, fp) ;    % get data

    [lightCommand_Exc(a,:), error] = ITCReadEpochStm(epochsExc(a), 0,fp);
    
    [SIexc, error] = ITCGetSamplingInterval(epochsExc(a), fp); % get sampling interval
    
    [seedExc(a),error] = ITCGetSeed(epochsExc(a),0,0,fp) ;
end

for a = 1:length(epochsInh) ; % for each g epoch
    [dataInh(a,:),error] = ITCReadEpoch(epochsInh(a), 0, fp) ;

    [lightCommand_Inh(a,:), error] = ITCReadEpochStm(epochsInh(a), 0,fp);
 
    [SIinh, error] = ITCGetSamplingInterval(epochsInh(a), fp);
    
    [seedInh(a),error] = ITCGetSeed(epochsInh(a),0,0,fp) ;
end


% for a = 1:length(epochsCA) ; % for each spike epoch
%     [dataCA(a,:),error] = ITCReadEpoch(epochsCA(a), 0, fp) ;    % get data
% 
%     [lightCommand_CA(a,:), error] = ITCReadEpochStm(epochsCA(a), 0,fp);
%     
%     [SIca, error] = ITCGetSamplingInterval(epochsCA(a), fp); % get sampling interval
%     
%     [seedCA(a),error] = ITCGetSeed(epochsCA(a),0,0,fp) ;
% end


% check sample intervals match
if SIexc~=SIinh && SIca ;
    display('exc and inh SI do not match')
else
    SI = SIexc * 1e-6; % Sampling interval in sec
    if Input(A).ITC18flag == 1 ;
        SI = SI*1.25 ;
    end
end

[prePnts, error] = ITCGetStmPrePts(epochsExc(1), 0, 0, fp) ; % points collected beyond which data is not worth analyzing
[postPnts, error] = ITCGetStmTailPts(epochsExc(1), 0, 0, fp) ;
[stmPnts, error] = ITCGetStmPts(epochsExc(1),0,0, fp) ;

% [prePntsCA, error] = ITCGetStmPrePts(epochsCA(1), 0, 0, fp) ; % points collected beyond which data is not worth analyzing
% [postPntsCA, error] = ITCGetStmTailPts(epochsCA(1), 0, 0, fp) ;
% [stmPntsCA, error] = ITCGetStmPts(epochsCA(1),0,0, fp) ;

samplerate = 1/SI ; % Hz at which data was collected
time = [SI:SI:SI*length(dataExc(1,:))] ; % time vector in seconds

% low pass filter to remove electrical crap
dataExc_lpf = lowPassFilter(dataExc, samplerate, 5000) ; %(signal,samplerate,cutoff frequ (hz))
dataInh_lpf = lowPassFilter(dataInh, samplerate, 5000) ; 

% high pass filter to remove slow drift 
% greg S. sent this to me to help implement a butterworth and avoid ringing
F=1 ; % filter cuttoff
Wn = F*SI(1); %normalized frequency cutoff
[z, p, k] = butter(1,Wn,'high'); %
[sos,g]=zp2sos(z,p,k); 
myfilt=dfilt.df2sos(sos,g);

dataExc_hpf = filter(myfilt,dataExc_lpf')'; % filter implementation
dataInh_hpf = filter(myfilt,dataInh_lpf')'; 

% change to conductances and subtract off means 
for a = 1:size(dataExc,1) ; % for each trial
    GExc_hpf(a,:) = dataExc_hpf(a,:)/-61 - mean(dataExc_hpf(a,prePnts:end-postPnts)/-61) ; % get conductance from stable currrents
end
for a = 1:size(dataInh,1) ; % for each trial
    GInh_hpf(a,:) = dataInh_hpf(a,:)/61 - mean(dataInh_hpf(a,prePnts:end-postPnts)/61) ; 
end

offsetExc = min(min(GExc_hpf(:,prePnts:end-postPnts))) ;
offsetInh = min(min(GInh_hpf(:,prePnts:end-postPnts))) ;

% ofset G (assumes all g have same mean and 1 min) 
for a = 1:size(dataExc,1) ; % for each trial
    GExc_hpf(a,:) = GExc_hpf(a,:) - offsetExc ; % offsets
end
for a = 1:size(dataInh,1) ;
    GInh_hpf(a,:) = GInh_hpf(a,:) - offsetInh ; 
end

% averages for each unique stimulus
UniqueSeeds = unique([seedExc,seedInh]) ;

for a=1:length(UniqueSeeds) ;   
    uiExc = find(seedExc==UniqueSeeds(a)) ;
    uiInh = find(seedInh==UniqueSeeds(a)) ;
    
    GExc{a} = GExc_hpf(uiExc,:) ;
    GInh{a} = GInh_hpf(uiInh,:) ;
    
    GExc_mean(a,:) = mean(GExc{a},1) ;
    GInh_mean(a,:) = mean(GInh{a},1) ;
    
    ccG(a,:) = xcov(GExc_mean(a,prePnts:end-postPnts),GInh_mean(a,prePnts:end-postPnts),'coeff') ; % cross correlation of mean g

end

ccG_mean = mean(ccG) ;
ccTime = ([1:length(ccG)] - (length(ccG)+1)/2)*SI ;


% figures
figure
plot(ccTime,ccG_mean)

% for Igor

identifier = ['cc',num2str(A)] ;
ForIgor.(identifier) = ccG_mean ;

identifier = ['ccTime',num2str(A)] ;
ForIgor.(identifier) = ccTime ;





