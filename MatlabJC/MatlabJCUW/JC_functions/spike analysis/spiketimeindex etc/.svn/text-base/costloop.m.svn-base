cd ~/data_analysis/Index ;

% Get appropriate cell structure format
CellInfo = LoadSCIData(CellInfo,1) ;  % for two amps used 

% this function excludes all the epochs but the ones you want
% Input = CellInfo, EpochConditionNumber, vector of the epochs you want
% Output = CellInfo, b = index of epochs in cellinfo
CellInfo = EpochExcluder(CellInfo,epochCond_num,epochsUwant) ; 

% this must come after EpochExcluder
EpochCondition = LoadAndSmoothEpochCondition(CellInfo,5000,1) ;

[SpikeTimeIndex]=GetSpikeTimes_CellAttachedjc(CellInfo,EpochCondition(epochCond_num), pre, post, threshold) ;

a=0 ;
for qt = 1:100:1001
    cost = 2/qt ;
    a=a+1
% This is the spike distance metric function as assessed by a cumulative
% probability histogram
[DeltaT, tli_spike, tlj_spike, Percent_Pairs_Quantified, DeltaT_Histogram,X_Values,DeltaT_CumProb] = DeltaT_Distribution_From_SCR(SpikeTimeIndex,cost) ;

% find median value of ditribution
[c,d] = min(abs(DeltaT_CumProb-.5)) ; % c=min value, d=index of min value
median_DeltaT(a) = X_Values(d) ; % what detlaT value is that? 

end

figure, plot([1:100:1001],median_DeltaT) ;
