%--------------------------------------------------------------------------
% TODO:
% - cell attached data for midget-parasol pairs
% - summarize cone noise measurements, fit to estimate noise correlation
% - summarize cone-rgc pairs
% - try models with (1) hexagonal cone lattice; (2) gaussian weights
% - check sizes of diffuse and midget bipolar axon terminals
%--------------------------------------------------------------------------

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

% modulated + constant light, on parasol pairs
load 071007Ac2      % coupled, very nice, figure
load 071007Ac3      % coupled, nice spike distance fig)
load 021307Hc3      % low/high light, nice spike distance at low/high light, nice mod
load 062707Ac2      % very nice, coupled, low overlap image fig, victor fig
load 091807Ac1      % possible figure for single cell response props

% modulated + constant light, off parasol pairs
load 071707Ac2      
load 080107Ac1
load 080107Ac2      % very nice - figure for vc, coupling, victor
load 080107Ac3      % very nice 
load 101007Ac1      % very nice - potential figure
load 101007Ac2

% direct measures of coupling
load 062707Ac2      % on parasol pair - very nice, coupled
load 070307Ac1      % on parasol pair - cc/vc
load 071007Ac2      % on pair, mod + contrast, coupled, very nice, figure
load 071007Ac3      % on pair, mod + contrast, coupled
load 091807Ac1

load 080107Ac2      % off pair - lots of dendritic overlap - images ok (check channel 2)
load 080107Ac3      % off pair
load 080107Ac1      % off pair
load 071707Ac2      % off pair
load 071807Ac2      % off pair, lousy light responses
load 071807Ac3      % off pair
load 072407Ac3      % off pair
load 101007Ac1      % very nice - potential figure
load 101007Ac2

% current clamp
load 092706Hc1      % cc/on cell pair, nice
load 092706Hc2      % cc/on cell pair, nice
load 110806Hc1      % cc/on cell pair, nice, cell attached; cell attached/cc 
load 110806Hc2      % cc/on cell pair, nice, cell attached; cell attached/cc; cc/cc

% on cell only
load 091306Hc1      % on cell only, nice image
load 092006Hc2      % on cell only
load 092006Hc3      % on cell only
load 092006Hc4      % on cell only
load 062707Ac1      % on cell only, mod + constant
load 071007Ac1      % mod + contrast, on pair
load 071807Ac1      % mod + contrast, off pair, step figure cell?
load 072407Ac2      % mod + contrast, on pair


%%

%--------------------------------------------------------------------------
% data loading
%--------------------------------------------------------------------------

%%
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

% plot average of each EpochCondition; threshold if on cell
Threshold = 0.2;
Polarity = -1;
DecimatePts = 4;
figure(1); clf
for cond = 1:length(EpochCondition)
    if (strcmp(EpochCondition(cond).Label, 'on cell'))
        clear PreData;
        PreSpikeTimes = GetSpikeTimes_CellAttached(CellInfo, EpochCondition(cond), Threshold, Polarity);
        PreData(1:size(EpochCondition(cond).EpochData.Data, 2)) = 0;
        for epoch = 1:size(PreSpikeTimes, 2)
            PreData(PreSpikeTimes{epoch} * 10) = PreData(PreSpikeTimes{epoch} * 10) + 1;
        end
        PreData = PreData / size(PreSpikeTimes, 2); 
        SmoothedPSTH = decimate(PreData, DecimatePts);
        plot(PreData);
    else
        plot(EpochCondition(cond).AverageResponse);
    end
    pause(1)
end

%%

%--------------------------------------------------------------------------
% correlations between spike times
%--------------------------------------------------------------------------

%%
label = 'on cell';
Threshold = 0.25;
PrePolarity = 1;
PostPolarity = 1;
DecimatePts = 4;
verbose = 1;
MaxLag = 2000;
PreOnCellFlag = 1;
PostOnCellFlag = 1;

clear FoundCond;
for cond = 1:length(EpochCondition)
    if (strcmp(EpochCondition(cond).Label, label))
        if (EpochCondition(cond).SegNum == 0)
            if (exist('FoundCond'))
                FoundCond = [FoundCond cond];
            else
                FoundCond = cond;
            end
        end
    end
end

PrePts = FindSearchPara(EpochCondition(FoundCond(1)), 'PrePoints');
StmPts = FindSearchPara(EpochCondition(FoundCond(1)), 'StimDur');
EndPts = length(EpochCondition(FoundCond(1)).AverageResponse);
StartPts = PrePts + StmPts + 5000;

if (length(FoundCond) > 0)
    for cond = 1:length(FoundCond)
    
        PreSpikeTimes = GetSpikeTimes(CellInfo, EpochCondition(FoundCond(cond)), Threshold, StartPts, EndPts, PreOnCellFlag, PrePolarity);
        PostSpikeTimes = GetSpikeTimes(CellInfo, EpochCondition(FoundCond(cond)+1), Threshold, StartPts, EndPts, PostOnCellFlag, PostPolarity);
        EndPts = size(EpochCondition(FoundCond(cond)).EpochData.Data, 2);

        clear CrossCorrelation;
        for epoch = 1:size(PreSpikeTimes, 2)

            PreData(1:size(EpochCondition(FoundCond(cond)).EpochData.Data, 2)) = 0;
            PostData(1:size(EpochCondition(FoundCond(cond)).EpochData.Data, 2)) = 0;
            PreData(PreSpikeTimes{epoch} * 10) = 1;
            PostData(PostSpikeTimes{epoch} * 10) = 1;
            
            if (verbose)
                fprintf(1, '%d %d %d\n', epoch, sum(PreData), sum(PostData));
            end

            if (mean(PreData) > 0)
                if (mean(PostData) > 0)             
                    if (exist('CrossCorrelation'))
                        CrossCorrelation = CrossCorrelation + xcorr(PreData, PostData, MaxLag, 'coeff');
                        PreAutoCorrelation = PreAutoCorrelation + xcorr(PreData, PreData, MaxLag, 'coeff');
                        PostAutoCorrelation = PostAutoCorrelation + xcorr(PostData, PostData, MaxLag, 'coeff');
                    else
                        CrossCorrelation = xcorr(PreData, PostData, MaxLag, 'coeff');
                        PreAutoCorrelation = xcorr(PreData, PreData, MaxLag, 'coeff');
                        PostAutoCorrelation = xcorr(PostData, PostData, MaxLag, 'coeff');
                    end
                end
            end
            
        end

        CrossCorrelation = CrossCorrelation / size(PreSpikeTimes, 2);
        PreAutoCorrelation = PreAutoCorrelation / size(PreSpikeTimes, 2);
        PostAutoCorrelation = PostAutoCorrelation / size(PostSpikeTimes, 2);
        % scale to firing rate
        CrossCorrelation = CrossCorrelation / 0.0001;
        tme = 1:length(CrossCorrelation);
        tme = (tme - length(CrossCorrelation)/2) * 0.0001;
        % fit with Gaussian
        coef = [10 0 0.01 1];
        fitcoef = nlinfit(tme(length(CrossCorrelation)/2 - 100:length(CrossCorrelation)/2 + 100), CrossCorrelation(length(CrossCorrelation)/2 - 100:length(CrossCorrelation)/2 + 100), 'gaussianplusmean', coef);
        fit = gaussianplusmean(fitcoef, tme);
        fprintf(1, 'width = %d\n', fitcoef(3));
        OnCellWidth(cond) = fitcoef(3);
        SmoothedCorrelation = decimate(CrossCorrelation, DecimatePts);
        SmoothedTme = 1:length(SmoothedCorrelation);
        SmoothedTme = (SmoothedTme - length(SmoothedTme)/2) * 0.0001 * DecimatePts;
        figure(1); clf
        set(1, 'Position', [1 scrsz(4) scrsz(3)/3 scrsz(4)/2.5]);
        plot(SmoothedTme, SmoothedCorrelation, 'LineWidth', 2);
        xlim([-0.04 0.04])
        xlabel('time (sec)');
        ylabel('firing rate');
        FileName = strcat(strcat(BaseFileName, 'on-cell'), num2str(FoundCond(cond)));
        saveas(gcf, FileName, 'pdf');
        figure(2); clf
        set(2, 'Position', [scrsz(3)/3 scrsz(4) scrsz(3)/3 scrsz(4)/2.5]);
        plot(tme, PreAutoCorrelation, tme, PostAutoCorrelation, 'k');
        xlim([0 0.02]);
        ylim([0 0.01]);
        % coupling efficiency
        CouplingWidth = OnCellWidth(cond) * 2 / 0.0001;
        CorrArea = sum(CrossCorrelation(length(CrossCorrelation)/2-CouplingWidth:length(CrossCorrelation)/2+CouplingWidth));
        CorrArea = CorrArea - sum(CrossCorrelation(1:2*CouplingWidth));
        CorrArea = CorrArea * 0.0001;
        fprintf(1, 'coupling efficiency = %d\n', CorrArea);
        CellParameters.CellAttachedXCor(cond, :) = CrossCorrelation;
        CellParameters.CellAttachedPreACor(cond, :) = PreAutoCorrelation;
        CellParameters.CellAttachedPostACor(cond, :) = PostAutoCorrelation;
    end
end


%%

%--------------------------------------------------------------------------
% voltage/current clamp correlations 
%--------------------------------------------------------------------------

%%
label = 'vc';
clear FoundCond;
MaxLag = 2000;

PostThresholdFlag = 0;
PreThresholdFlag = 0;

for cond = 1:length(EpochCondition)
    if (strncmp(EpochCondition(cond).Label, label, length(label)))
        if (EpochCondition(cond).SegNum == 0)
            if (exist('FoundCond'))
                FoundCond = [FoundCond cond];
            else
                FoundCond = cond;
            end
        end
    end
end

PrePts = FindSearchPara(EpochCondition(FoundCond(1)), 'PrePoints');
StmPts = FindSearchPara(EpochCondition(FoundCond(1)), 'StimDur');
EndPts = length(EpochCondition(FoundCond(1)).AverageResponse);
StartPts = PrePts + StmPts + 5000;

if (length(FoundCond) > 0)
    for cond = 1:length(FoundCond)
    
        if (PreThresholdFlag)
            PreSpikeTimes = GetSpikeTimes(CellInfo, EpochCondition(FoundCond(cond)), Threshold);
        end
        if (PostThresholdFlag)
            PostSpikeTimes = GetSpikeTimes(CellInfo, EpochCondition(FoundCond(cond)+1), Threshold);
        end

        PreEpochData = GetGoodEpochData(CellInfo, EpochCondition(FoundCond(cond)));
        PostEpochData = GetGoodEpochData(CellInfo, EpochCondition(FoundCond(cond)+ 1));

        EndPts = size(PreEpochData, 2);

        clear CrossCorrelation;
        figure(2)
            
        for epoch = 1:size(PostEpochData, 1)

            if (PreThresholdFlag)
                PreData(1:size(EpochCondition(FoundCond(cond)).EpochData.Data, 2)) = 0;
                PreData(PreSpikeTimes{epoch} * 10) = 1;
                PreData = PreData(StartPts:EndPts);
            else
                PreData = PreEpochData(epoch, StartPts:EndPts) - mean(PreEpochData(:, StartPts:EndPts));
                PreData = PreData / max(abs(PreData));
                PreData = PreData - mean(PreData);
            end
            if (PostThresholdFlag)
                PostData(1:size(EpochCondition(FoundCond(cond)+1).EpochData.Data, 2)) = 0;
                PostData(PostSpikeTimes{epoch} * 10) = 1;
                PostData = PostData(StartPts:EndPts);
            else
                PostData = PostEpochData(epoch, StartPts:EndPts) - mean(PostEpochData(:, StartPts:EndPts));
                PostData = PostData / max(abs(PostData));
                PostData = PostData - mean(PostData);
            end

            fprintf(1,'variance = %d %d\n', var(PreData), var(PostData));
             
            if (exist('CrossCorrelation'))
                CrossCorrelation = CrossCorrelation + xcorr(PreData, PostData, MaxLag, 'coeff');
            else
                CrossCorrelation = xcorr(PreData, PostData, MaxLag, 'coeff');
            end
        end
        CrossCorrelation = CrossCorrelation / size(PostEpochData, 1);
        tme = 1:length(CrossCorrelation);
        tme = (tme - length(CrossCorrelation)/2) * 0.0001;
        figure(1); clf
        plot(tme, xcorr(mean(PreEpochData(:, PrePts:StmPts)), mean(PostEpochData(:, PrePts:StmPts)), MaxLag, 'coeff'), 'LineWidth', 2);
        axis tight
        % fit with Gaussian
        coef = [1 0 0.01 0.1];
        fitcoef = nlinfit(tme(length(CrossCorrelation)/2 - 100:length(CrossCorrelation)/2 + 100), CrossCorrelation(length(CrossCorrelation)/2 - 100:length(CrossCorrelation)/2 + 100), 'gaussianplusmean', coef);
        fit = gaussianplusmean(fitcoef, tme);
        fprintf(1, 'width = %d\n', fitcoef(3));
        WholeCellCCWidth(cond) = abs(fitcoef(3));
        figure(2); clf
        set(2, 'Position', [scrsz(3)/3 scrsz(4) scrsz(3)/3 scrsz(4)/2.5]);
%        plot(tme, CrossCorrelation, tme, fit, 'r');
        plot(tme, CrossCorrelation, 'LineWidth', 2);
        xlim([-0.1 0.1]);
        xlabel('time (sec)');
        ylabel('Correlation');
        FileName = strcat(strcat(BaseFileName, 'xcor'), num2str(FoundCond(cond)));
        saveas(gcf, FileName, 'pdf');

        clear PreAutoCorrelation PostAutoCorrelation;
        figure(3)
        for epoch = 1:size(PostEpochData, 1)

            PreData = PreEpochData(epoch, StartPts:EndPts) - mean(PreEpochData(:, StartPts:EndPts));
            PostData = PostEpochData(epoch, StartPts:EndPts) - mean(PostEpochData(:, StartPts:EndPts));
            PreData = PreData / max(abs(PreData));
            PostData = PostData / max(abs(PostData));
            PreData = PreData - mean(PreData);
            PostData = PostData - mean(PostData);

            if (exist('PreAutoCorrelation'))
                PreAutoCorrelation = PreAutoCorrelation + xcorr(PreData, PreData, MaxLag, 'coeff');
                PostAutoCorrelation = PostAutoCorrelation + xcorr(PostData, PostData, MaxLag, 'coeff');
            else
                PreAutoCorrelation = xcorr(PreData, PreData, MaxLag, 'coeff');
                PostAutoCorrelation = xcorr(PostData, PostData, MaxLag, 'coeff');
            end
        end
        PreAutoCorrelation = PreAutoCorrelation / size(PostEpochData, 1);
        PostAutoCorrelation = PostAutoCorrelation / size(PostEpochData, 1);
        tme = 1:length(PreAutoCorrelation);
        tme = (tme - length(PreAutoCorrelation)/2) * 0.0001;
        % fit with Gaussian
        coef = [1 0 0.01 0.1];
        fitcoef = nlinfit(tme(length(PreAutoCorrelation)/2 - 100:length(PreAutoCorrelation)/2 + 100), PreAutoCorrelation(length(PreAutoCorrelation)/2 - 100:length(PreAutoCorrelation)/2 + 100), 'gaussianplusmean', coef);
        fit = gaussianplusmean(fitcoef, tme);
        PreACWidth(cond) = abs(fitcoef(3));
        fitcoef = nlinfit(tme(length(PostAutoCorrelation)/2 - 100:length(PostAutoCorrelation)/2 + 100), PostAutoCorrelation(length(PostAutoCorrelation)/2 - 100:length(PostAutoCorrelation)/2 + 100), 'gaussianplusmean', coef);
        fit = gaussianplusmean(fitcoef, tme);
        PostACWidth(cond) = abs(fitcoef(3));
        fprintf(1, 'ac width = %d %d\n', PreACWidth(cond), PostACWidth(cond));
        figure(3)
        set(3, 'Position', [2*scrsz(3)/3 scrsz(4) scrsz(3)/3 scrsz(4)/2.5]);
        plot(tme, PreAutoCorrelation, tme, PostAutoCorrelation, 'k', tme, fit, 'r');
%        plot(tme, PostAutoCorrelation, 'k', 'LineWidth', 2);
        xlim([-0.1 0.1]);
        xlabel('time (sec)');
        ylabel('Correlation');
        FileName = strcat(strcat(BaseFileName, 'acor'), num2str(FoundCond(cond)));
        saveas(gcf, FileName, 'pdf');
        fprintf(1, 'condition %d (epoch = %d)\n', cond, EpochCondition(FoundCond(cond)).EpochNumbers(1));
        EpochCondition(FoundCond(cond)).Label
        pause(0.1)
        CellParameters.XCor(cond, :) = CrossCorrelation;
        CellParameters.ACor1(cond, :) = PreAutoCorrelation;
        CellParameters.ACor2(cond, :) = PostAutoCorrelation;
    end
end

%%

%--------------------------------------------------------------------------
% Batch Analysis - paper figs
%--------------------------------------------------------------------------

%%
% define plot color sequence, axis fonts
PlotColors = 'bgrkymcbgrkymcbgrkymcbgrkymc';
set(0, 'DefaultAxesFontName','Helvetica')
set(0, 'DefaultAxesFontSize', 14)
colormap([0 0 0])
scrsz = get(0, 'ScreenSize');

%% figure 1


cd '~/analysis/MATLAB/correlated-activity';
ExampleOnCell = '071007Ac2';
ExampleOffCell = '080107Ac2';
startpt = 56000;
endpt = 67000;

figure(1);clf;subplot(3, 2, 1);
tme = 1:endpt-startpt;
tme = tme * 0.0001;

AnalysisFileName = strcat(ExampleOnCell, 'analysis');
load(AnalysisFileName);
plot(tme, CellParameters.stm(startpt+1:endpt) * 5000 + 1000);
hold on

plot(tme, CellParameters.RawData(1, startpt+1:endpt), 'k');
plot(tme, CellParameters.RawData(2, startpt+1:endpt)-2500, 'r');
axis tight
xlabel('sec');
ylabel('pA');

subplot(3, 2, 2)
tme = 1:length(CellParameters.XCorMod);
tme = (tme - length(tme) / 2) * 0.0001;

plot(tme, CellParameters.XCorMod, 'k', tme, CellParameters.XCorResidual, 'r', tme, CellParameters.XCorConstant, 'b', 'LineWidth', 2);
xlim([-0.05 0.05]);
ylim([-0.3 1]);
ylabel('correlation');
xlabel('sec');

subplot(3, 2, 3);
AnalysisFileName = strcat(ExampleOffCell, 'analysis');
load(AnalysisFileName);
tme = 1:endpt-startpt;
tme = tme * 0.0001;

hold on
plot(tme, CellParameters.RawData(1, startpt+1:endpt), 'k');
plot(tme, CellParameters.RawData(2, startpt+1:endpt)-1000, 'r');
axis tight
xlabel('sec');
ylabel('pA');

subplot(3, 2, 4)
tme = 1:length(CellParameters.XCorMod);
tme = (tme - length(tme) / 2) * 0.0001;
ylabel('correlation');
xlabel('sec');

plot(tme, CellParameters.XCorMod, 'k', tme, CellParameters.XCorResidual, 'r', tme, CellParameters.XCorConstant, 'b', 'LineWidth', 2);
xlim([-0.05 0.05]);
ylim([-0.1 1]);
ylabel('correlation');
xlabel('sec');

clear CellName XCor XCorMax Overlap;
OverlapThreshold = 10;

CellName(1).name = '071007Ac2';
CellName(2).name = '071007Ac3';
CellName(3).name = '062707Ac2';
CellName(4).name = '021307Hc3';
CellName(5).name = '020807Hc1';
CellName(6).name = '111506Hc3';
CellName(7).name = '121206Hc3';     % check images
CellName(8).name = '072506Hc3';
CellName(9).name = '091306Hc2';    % recheck overlap - contamination
CellName(9).name = '100406Hc2';
CellName(10).name = '091807Ac1';
%CellName(11).name = '092006Hc1';       % work on images

for pair = 1:length(CellName);
    CellName(pair).name
    clear CellParameters;
    AnalysisFileName = strcat(CellName(pair).name, 'analysis');
    load(AnalysisFileName);
    figure(2)
    if (size(CellParameters.XCor, 1) > 1)
        for cnt = 1:size(CellParameters.XCor, 1)
            plot(CellParameters.XCor(cnt, :))
            pause()
        end
        xcorcnt = input('which one?');
    else
        xcorcnt = 1;
    end
    figure(3)
    plot(CellParameters.XCor(xcorcnt, :));
    pause(1);
    XCor(pair, :) = CellParameters.XCor(xcorcnt, :);
    XCorMax(pair) = max(CellParameters.XCor(xcorcnt, :));
    Indices = find(CellParameters.XCor(xcorcnt, :) > XCorMax(pair)/2);
    Halfwidth(pair) = (max(Indices) - min(Indices)) * 0.1;
    OverlapIndices = find(CellParameters.ImageParameters.NNDistances < OverlapThreshold);
    Overlap(pair) = length(OverlapIndices) / length(CellParameters.ImageParameters.NNDistances);
end
MeanXCor = mean(XCor);
SEMXCor = std(XCor) / sqrt(length(CellName));

tme = 1:length(XCor);
tme = (tme - length(tme) / 2) * 0.0001;

figure(1);subplot(3, 2, 5)
plot(Overlap, XCorMax, 'o');
xlim([0 max(Overlap)]);
ylim([0 max(XCorMax)]);
xlabel('overlap');
ylabel('peak correlation');

clear CellName XCorMod XCorResidual XCorConstant XCorConstantMax XCorResidualMax;
CellName(1).name = '071007Ac2';
CellName(2).name = '071007Ac3';
CellName(3).name = '021307Hc3';
CellName(4).name = '062707Ac2';
CellName(5).name = '091807Ac1';

for pair = 1:length(CellName);
    AnalysisFileName = strcat(CellName(pair).name, 'analysis');
    load(AnalysisFileName);
    XCorConstantMax(pair) = max(CellParameters.XCorConstant);
    XCorResidualMax(pair) = max(CellParameters.XCorResidual);
end

subplot(3, 2, 6)
plot(XCorConstantMax, XCorResidualMax, 'o');
hold on


% off cell pairs
CellName(1).name = '071707Ac2';
CellName(2).name = '080107Ac2';
CellName(3).name = '080107Ac3';
CellName(4).name = '080107Ac1';
CellName(5).name = '101007Ac1';
CellName(6).name = '101007Ac2';

for pair = 1:length(CellName);
    AnalysisFileName = strcat(CellName(pair).name, 'analysis');
    load(AnalysisFileName);
    XCorConstantMax(pair) = max(CellParameters.XCorConstant);
    XCorResidualMax(pair) = max(CellParameters.XCorResidual);
    Indices = find(CellParameters.XCorConstant(xcorcnt, :) > XCorConstantMax(pair)/2);
    Halfwidth(pair) = (max(Indices) - min(Indices)) * 0.1;
end

subplot(3, 2, 6)
plot(XCorConstantMax, XCorResidualMax, 'or');
plot([0 1], [0 1]);
hold off
xlim([0 0.5]);
ylim([0 0.5]);
xlabel('peak constant');
ylabel('peak modulated');


%%
% tonic current during constant light
%

FreqCutoff = 2000;
label = 'vc';
DecimatePts = 1;
OnFlag = 0;

clear CellName
%  modulated + constant light, on parasol pairs
if (OnFlag)
    CellName(1).name = '071007Ac3';
    CellName(2).name = '021307Hc3';
    CellName(3).name = '062707Ac2';
    CellName(4).name = '091807Ac1';
    CellName(5).name = '021308Ac1';
    CellName(6).name = '021308Ac3';
    CellName(7).name = '032508Ac3';    
    CellName(8).name = '032708Ac2';
else
    % off cell pairs
    CellName(1).name = '071707Ac2';
    CellName(2).name = '080107Ac2';
    CellName(3).name = '080107Ac3';
    CellName(4).name = '080107Ac1';
    CellName(5).name = '101007Ac1';
    CellName(6).name = '101007Ac2';
    CellName(7).name = '032508Ac2';
    CellName(8).name = '040208Ac2';
    CellName(9).name = '040208Ac3';    
end

for cell = 1:length(CellName)
    
    IndexFolder; cd monk-rgc;
    load(CellName(cell).name);

    % read in and smooth EpochCondition data
    if (isfield(CellInfo, 'EpochCondition')) 
        CellInfo = LoadSCIData(CellInfo, 1);
        EpochCondition = LoadAndSmoothEpochCondition(CellInfo, FreqCutoff, DecimatePts);
    end

    AnalysisFolder; cd correlated-activity;
    clear CellParameters
    DataFileName = strcat(CellName(cell).name, 'analysis.mat');
    if (exist(DataFileName))
        load(DataFileName);
    end

    background = 1;
    clear FoundCond;

    for cond = 1:length(EpochCondition)
        if (strncmp(EpochCondition(cond).Label, label, length(label)))
            if (EpochCondition(cond).SegNum == 0)
                if (exist('FoundCond'))
                    FoundCond = [FoundCond cond];
                else
                    FoundCond = cond;
                end
            end
        end
    end

    for cond = 1:length(FoundCond)
        CellParameters.PreTonicOffset(cond) = max(EpochCondition(FoundCond(cond)).AverageResponse(40000:50000)) / -min(EpochCondition(FoundCond(cond)).AverageResponse(40000:50000));
        CellParameters.PostTonicOffset(cond) = max(EpochCondition(FoundCond(cond)+1).AverageResponse(40000:50000)) / -min(EpochCondition(FoundCond(cond)+1).AverageResponse(40000:50000));
    end
    
    save(DataFileName, 'CellParameters');
end

%%
OnFlag = 0;

clear CellName TonicOffset
if (OnFlag)
    CellName(1).name = '071007Ac3';
    CellName(2).name = '021307Hc3';
    CellName(3).name = '062707Ac2';
    CellName(4).name = '091807Ac1';
    CellName(5).name = '021308Ac1';
    CellName(6).name = '021308Ac3';
    CellName(7).name = '032508Ac3';    
    CellName(8).name = '032708Ac2';
else
    % off cell pairs
    CellName(1).name = '071707Ac2';
    CellName(2).name = '080107Ac2';
    CellName(3).name = '080107Ac3';
    CellName(4).name = '080107Ac1';
    CellName(5).name = '101007Ac1';
    CellName(6).name = '101007Ac2';
    CellName(7).name = '032508Ac2';
    CellName(8).name = '040208Ac2';
    CellName(9).name = '040208Ac3';    
end

for cell = 1:length(CellName)
     
    AnalysisFolder; cd correlated-activity;
    clear CellParameters
    DataFileName = strcat(CellName(cell).name, 'analysis.mat');
    if (exist(DataFileName))
        load(DataFileName);
    end
    
    TonicOffset(2*cell-1) = CellParameters.PreTonicOffset(CellParameters.xcorcnt);
    TonicOffset(2*cell) = CellParameters.PostTonicOffset(CellParameters.xcorcnt);
    
end

