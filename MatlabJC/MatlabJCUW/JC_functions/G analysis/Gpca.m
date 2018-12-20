
perform.PCAinh = 1 ;


for A=2 ;
%% 
Input(2).cellname = '060908Bc4' ; % 
Input(2).CA = '[45:74]' ;
Input(2).Exc = '[110:139]' ; % VC at inh reversal pot
Input(2).Inh = '[163:182]' ;
Input(2).InhApb = '[212:241]' ; % VC at exc revesal with APB present
Input(2).ExcApb = '[270:289]' ;
Input(2).InhWash = '[541:570]' ; % inh post apb
Input(2).ExcWash = '[616:645]' ;

Input(4).cellname = '061008Bc4' ; % 
Input(4).CA = '[42:81]' ;
Input(4).Exc = '[106:125]' ; % VC at inh reversal pot
Input(4).Inh = '[147:166]' ;
Input(4).InhApb = '[212:251]' ; % VC at exc revesal with APB present
Input(4).ExcApb = '[292:325]' ;

% Parameters assuming sample rate at 10 kHz
Parameters.PrePnts = 10000 ;    % number of points before stimuli starts
Parameters.StmPnts = 60000 ;    % number of points during stimuli
Parameters.PostPnts = 1000 ;    % number of points after stimuli ends
Parameters.STAPnts = 3000 ;     % number of points used for prespike wave forms
Parameters.DecimatePnts = 10 ;  % number of points averaged together in order to downsample prespike waveforms
Parameters.SmoothPnts = 100 ;   % number of points used to smooth spike train to make PSTH
Parameters.ChopPnts = 6000 ;

%%
[fp, error] = ITCInitializeAnalysis(500000,  ['~/Data/Primate/',Input(A).cellname]);

%% 

epochs = str2num(Input(A).Exc) ;
for a = 1:length(epochs) ; % for each spike epoch
    [stm(a,:), error] = ITCReadEpochStm(epochs(a), 0, fp);  % get the light stimulus
    [data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp);    % get cell attached data
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp);
    SR(a) = 1/(SI(a) * 1e-6); % Sampling rate in Hz
   
end

[prePnts, error] = ITCGetStmPrePts(epochs(1), 0, 0, fp) ; % points collected beyond which data is not worth analyzing
[postPnts, error] = ITCGetStmTailPts(epochs(1), 0, 0, fp) ;

data = data - repmat(mean(data,2),1,size(data,2)) ; % subtract mean 

% rectifier, zero, and normalize stim to a variance of 1
negstmPnts = find(stm<0) ;  % find indicies which would be getting a negative stim voltage
stm(negstmPnts) = 0 ;       % make these points zero
stm = stm - repmat(mean(stm(:,prePnts:end-postPnts),2),1,size(stm,2)) ; % subtract off mean of stimulus during time varying stimulus
stm = stm./repmat(std(stm(:,Parameters.PrePnts:Parameters.PrePnts+Parameters.StmPnts),[],2),1,size(stm,2)) ;


% covariance
MaxLag = Parameters.STAPnts/Parameters.DecimatePnts ;
PrePts = Parameters.PrePnts ;
StmPts = Parameters.StmPnts  ;

StmCovar(1:MaxLag, 1:MaxLag) = 0;
RespCovar(1:MaxLag, 1:MaxLag) = 0;

for trial = 1:size(data,1) ;
        ModData(trial,:) = DecimateWave(data(trial,PrePts:PrePts + StmPts), Parameters.DecimatePnts);
        ModStm(trial,:) = DecimateWave(stm(trial,PrePts:PrePts + StmPts), Parameters.DecimatePnts);
end

[covar,EigVal,EigVec] = contPCA(ModData,ModStm,MaxLag) ;

filter1 = fliplr(EigVec(:,1)') ;
filter1 = interp1(filter1,[1:.1:300]) ;

filter2 = fliplr(EigVec(:,2)') ;
filter2 = interp1(filter2,[1:.1:300]) ;

for trial = 1:size(data,1) ;
    LinPredict1(trial,:) = conv(stm(trial,:),filter1) ;
    LinPredict2(trial,:) = conv(stm(trial,:),filter2) ;
end

[InputBins1,InputBins2,Nonlinearity,Nonlinearity_sem] = NonLinFilterFinder2(LinPredict1(:,PrePts:PrePts + StmPts),LinPredict2(:,PrePts:PrePts + StmPts),data(:,PrePts:PrePts + StmPts),10,10) ;
Predict = NonLinFilter2(InputBins2,InputBins1,LinPredict2(1,PrePts:PrePts + StmPts),LinPredict1(1,PrePts:PrePts + StmPts),Nonlinearity) ;


figure
plot(LinPredict(1,PrePts:PrePts + StmPts),data(1,PrePts:PrePts + StmPts))
[InputBins,Nonlinearity,Nonlinearity_sem] = NonLinFilterFinder(LinPredict(1,PrePts:PrePts + StmPts),data(1,PrePts:PrePts + StmPts),10) ;
hold on
plot(InputBins,Nonlinearity,'r')
plot(InputBins,Nonlinearity+Nonlinearity_sem,'r--')
plot(InputBins,Nonlinearity-Nonlinearity_sem,'r--')

Predict = NonLinFilter2(LinPredict2,InputBins,Nonlinearity) ;


epochs = str2num(Input(A).InhApb) ;
for a = 1:length(epochs) ; % for each spike epoch
    [stm(a,:), error] = ITCReadEpochStm(epochs(a), 0, fp);  % get the light stimulus
    [data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp);    % get cell attached data
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp);
    SR(a) = 1/(SI(a) * 1e-6); % Sampling rate in Hz
end

dataAPB = data - repmat(mean(data,2),1,size(data,2)) ; % subtract mean 

end
