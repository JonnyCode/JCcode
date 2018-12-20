% this script will create normalized tuning curve for export to igor for
% suplemental figure comparing cell attached data and dynamic clamp data
% 7/20/10

A=85 ;
id = 'CADS' ;
epochs = str2num(Input(A).(id)) ;

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
numTrials = length(epochs) ;

for a = 1:numTrials ; % for each spike epoch
    [dataCA1(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data
   
    [SI_CA(a), error] = ITCGetSamplingInterval(epochs(a), fp); % get sampling interval
    SI_CA(a) = SI_CA(a) * 1e-6; % Sampling interval in sec
end

if Input(A).ITC18flag == 1 ;
    SI_CA = SI_CA*1.25 ;
end

time_CA = [1:length(dataCA1)]*SI_CA(1) ;

SpatialStimParams = hdf5load(['~/Data/mouse/',Input(A).cellname,'_spatial.h5']) ; % load spatial stim params
frameRate = 60 ;

% spike detections
SpikePnts1 = SpikeDetection(dataCA1,10,1/SI_CA(1)) ;

SpikeTrain1 = zeros(size(dataCA1)) ;

for a=1:numTrials ;
    SpikeTrain1(a,SpikePnts1{a}) = 1 ;
end

% bar spikes
for a = 1:numTrials ;
    StrucString = ['params_epoch',num2str(epochs(a))] ; 
    Struct = SpatialStimParams.(StrucString) ;

    BarAngles(a,:) = Struct.BarAngle ;
    PrePts = floor(Struct.spatial_prepts/(60*SI_CA(a))) ;
    PostPts = floor(Struct.spatial_postpts/(60*SI_CA(a))) ;

    GroupPts = floor((length(SpikeTrain1)-PrePts-PostPts)/size(BarAngles,2)) ;

    spikeRate(a,:) = DecimateWave(SpikeTrain1(a,PrePts:end-PostPts), GroupPts)./ SI_CA(a) ;

    [SortedBarAngles,i] = sort(BarAngles(a,:)) ;
    SortedSpikeRate(a,:) = spikeRate(a,i) ;

end

FR_mean = mean(SortedSpikeRate) ;

SortedSpikeRate_norm = SortedSpikeRate/max(FR_mean) ;

FR_mean_norm = mean(SortedSpikeRate_norm) ;
FR_std_norm = std(SortedSpikeRate_norm) ;


% FORIGOR
identifier = ['meanFRCANormLin',num2str(A)] ; % added 7/20/10
ForIgor.(identifier) = FR_mean_norm ;

identifier = ['stdFRCANormLin',num2str(A)] ; % added 7/20/10
ForIgor.(identifier) = FR_std_norm ;



