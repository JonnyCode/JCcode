%%
% define plot color sequence, axis fonts

PlotColors = 'bgrkymcbgrkymcbgrkymcbgrkymc';
set(0, 'DefaultAxesFontName','Helvetica')
set(0, 'DefaultAxesFontSize', 22)
colormap([0 0 0])
scrsz = get(0, 'ScreenSize');

%%
%--------------------------------------------------------------------------
% data sets
IndexFolder; cd 'monk-rgc';

%%
%%--------------------------------------------------------------------------
% data loading
%--------------------------------------------------------------------------

cd '~/analysis/MATLAB/correlated-activity';

% check if local or server account
DataFolderPath = '~/Data/';
DataFileName = CellInfo.CellFile;
Indices = strfind(DataFileName, '/');
BaseFileName = DataFileName(max(Indices)+1:length(DataFileName));
CellInfo.CellFile = strcat(DataFolderPath, BaseFileName);

clear CellParameters
DataFileName = strcat(BaseFileName, 'analysis.mat');
if (exist(DataFileName))
    load(DataFileName);
end

% read in and smooth EpochCondition data
FreqCutoff = 2000;
DecimatePts = 1;
if (isfield(CellInfo, 'EpochCondition')) 
    CellInfo = LoadSCIData(CellInfo, 1);
    EpochCondition = LoadAndSmoothEpochCondition(CellInfo, FreqCutoff, DecimatePts);
end
if (isfield(CellInfo, 'FamilyCondition'))
    CellInfo = LoadSCIData(CellInfo, 2);
    FreqCutoff = 2000;
    FamilyCondition = LoadAndSmoothFamilyCondition(CellInfo, FreqCutoff, DecimatePts);
end

%%
%----------------------------------------------------------------------------------
% Modulation from CellInfo
%----------------------------------------------------------------------------------

[fp, error] = ITCInitializeAnalysis(500000, CellInfo.CellFile);

%% 
% covariance
MaxLag = 100;
cond = 1;
DecimatePts = 40;

PrePts = FindSearchPara(EpochCondition(cond), 'PrePoints');
StmPts = FindSearchPara(EpochCondition(cond), 'StimDur');
SamplingInterval = FindSearchPara(EpochCondition(cond), 'SampInterv') * 1e-6;
NumEpochs = 0;
clear StmCovar RespCovar;
StmCovar(1:MaxLag, 1:MaxLag) = 0;
RespCovar(1:MaxLag, 1:MaxLag) = 0;

for epoch = 1:length(EpochCondition(cond).EpochNumbers)
    if (EpochCondition(cond).ExcludeEpochs(epoch) ~= 1)
        [data, error] = ITCReadEpoch(EpochCondition(cond).EpochNumbers(epoch), 0, fp);
        [stm,error] = ITCReadEpochStm(EpochCondition(cond).EpochNumbers(epoch), 0, fp);
        ModData = decimate(data(PrePts:PrePts + StmPts), DecimatePts);
        ModStm = decimate(stm(PrePts:PrePts + StmPts), DecimatePts);
        ModStm = ModStm - mean(ModStm);
        ModStm = ModStm / var(ModStm);
        ModData = ModData - mean(ModData);
        ITCAnalogCovar(length(ModStm), MaxLag, ModStm, ModData, RespCovar, StmCovar);
    end
end

% symmetrize covariance matrix
for i = 1:MaxLag
    for j = 1:i
        RespCovar(i, j) = RespCovar(j, i);
    end
end

figure(1);
mesh(RespCovar);
%%
% PCA

[EigVec, EigVal] = eig(RespCovar);
figure(2);
semilogy(abs(EigVal), 'o');

[EigVec, EigVal] = eigs(RespCovar, 2);

figure(3);
plot(EigVec);

%%

