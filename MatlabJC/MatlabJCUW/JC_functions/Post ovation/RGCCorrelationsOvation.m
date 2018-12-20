%%
% general setup stuff

% define plot color sequence, axis fonts
PlotColors = 'bgrkymcbgrkymcbgrkymcbgrkymc';
set(0, 'DefaultAxesFontName','Helvetica')
set(0, 'DefaultAxesFontSize', 18)
colormap([0 0 0])
scrsz = get(0, 'ScreenSize');
%%
% database communication
import auimodel.*

%%
% on midget pairs
cd ~/Data/Ovation-Data;
load on-midget-pairs.mat
epochs = EpochList(OvationExport);
cd ~/Analysis/MATLAB/correlated-activity;

%%
% on/off parasol pairs
cd ~/Data/Ovation-Data;
load on-off-parasol-pairs.mat
epochs = EpochList(OvationExport);
cd ~/Analysis/MATLAB/correlated-activity;

%%
% create EpochTree
tree = EpochTree(epochs, {'protocolSettings.acquirinoCellBasename', 'protocolSettings.PreSynapticHold', 'protocolSettings.PostSynapticHold'});
visualize(tree);

%%
% look at average responses for each leaf of tree
figure(1); clf;
for leaf = 1:length(tree.leafNodes)
    PreEpochData = tree.leafNodes{leaf}.epochList.responsesByStreamName('Amp_1');
    PostEpochData = tree.leafNodes{leaf}.epochList.responsesByStreamName('Amp_2');
    plot([1:size(PreEpochData, 2)], mean(PreEpochData), [1:size(PostEpochData, 2)], mean(PostEpochData));
    pause(1);
end 

%%
% crosscorrelation functions
MaxLag = 2000;                  % time window for cross-correlation function
ResponseDelayPts = 5000;        % exclude at and of modulated light    
ExcHold = -55;                  % typical excitatory holding potential
InhHold = +25;                  % typical inhibitory holding potential

% find holding potentials for each leaf
for leaf = 1:length(tree.leafNodes)
    SampleEpoch = tree.leafNodes{leaf}.epochList.elements{1};
    HoldString = SampleEpoch.protocolSettings.PreSynapticHold;
    Index = find(HoldString == 'm');
    PreHold(leaf) = str2num(HoldString(1:Index-1));
    HoldString = SampleEpoch.protocolSettings.PostSynapticHold;
    Index = find(HoldString == 'm');
    PostHold(leaf) = str2num(HoldString(1:Index-1));
end

% find exc-exc conditions
Leaves{1} = find(abs(PreHold - ExcHold) < 20 & abs(PostHold - ExcHold) < 20);

% find inh-inh conditions
Leaves{2} = find(abs(PreHold - InhHold) < 20 & abs(PostHold - InhHold) < 20);

% find exc-inh and inh-exc conditions
tempLeaves = find(abs(PreHold - ExcHold) < 20 & abs(PostHold - InhHold) < 20);
Leaves{3} = [tempLeaves find(abs(PreHold - InhHold) < 20 & abs(PostHold - ExcHold) < 20)];

% loop across each set of conductances as defined above
for cond = 1:length(Leaves)
    
    % loop across each lead in given set
    for leaf = 1:length(Leaves{cond})

        % use sample epoch to extract key parameters for epochs in leaf
        SampleEpoch = tree.leafNodes{Leaves{cond}(leaf)}.epochList.elements{1};
        StimulusString = strcat('stimuli:', tree.leafNodes{Leaves{cond}(leaf)}.epochList.stimuliStreamNames{1});
        fprintf(1, 'starting %s PreHold %s PostHold %s\n', SampleEpoch.protocolSettings.acquirinoCellBasename, SampleEpoch.protocolSettings.PreSynapticHold, SampleEpoch.protocolSettings.PostSynapticHold);
        PrePts = SampleEpoch.protocolSettings.(strcat(StimulusString,':prepts'));
        StmPts = SampleEpoch.protocolSettings.(strcat(StimulusString,':stmpts'));
        EndPts = SampleEpoch.protocolSettings.(strcat(StimulusString, ':tailpts'));
        SamplingInterval = SampleEpoch.protocolSettings.('acquirinoSamplingInterval') * 1e-6;

        % start of constant light analysis
        StartPts = PrePts + StmPts + ResponseDelayPts;

        % get data
        PreEpochData = tree.leafNodes{Leaves{cond}(leaf)}.epochList.responsesByStreamName('Amp_1');
        PostEpochData = tree.leafNodes{Leaves{cond}(leaf)}.epochList.responsesByStreamName('Amp_2');

        EndPts = size(PreEpochData, 2);

        clear CrossCorrelation;

        % loop across all epochs
        for epoch = 1:size(PostEpochData, 1)

            PreData = PreEpochData(epoch, StartPts:EndPts) - mean(PreEpochData(:, StartPts:EndPts));
            PreData = PreData / max(abs(PreData));
            PreData = PreData - mean(PreData);

            PostData = PostEpochData(epoch, StartPts:EndPts) - mean(PostEpochData(:, StartPts:EndPts));
            PostData = PostData / max(abs(PostData));
            PostData = PostData - mean(PostData);

            if (exist('CrossCorrelation'))
                CrossCorrelation = CrossCorrelation + xcorr(PreData, PostData, MaxLag, 'coeff');
            else
                CrossCorrelation = xcorr(PreData, PostData, MaxLag, 'coeff');
            end
        end
        CrossCorrelation = CrossCorrelation / size(PostEpochData, 1);
        tme = 1:length(CrossCorrelation);
        tme = (tme - length(CrossCorrelation)/2) * SamplingInterval;
        
        % plot and save images
        figure(1); clf
        set(1, 'Position', [1 scrsz(4) scrsz(3)/3 scrsz(4)/2.5]);
        PreMean = mean(PreEpochData(:, PrePts:StmPts));
        PostMean = mean(PostEpochData(:, PrePts:StmPts));
        plot([1:length(PreMean)], PreMean, [1:length(PostMean)], PostMean);
        xlim([10000 25000]);
        title(SampleEpoch.protocolSettings.acquirinoCellBasename);
        legend(SampleEpoch.protocolSettings.PreSynapticHold, SampleEpoch.protocolSettings.PostSynapticHold, 'Location', 'Best');
        legend boxoff

        figure(2); clf
        set(2, 'Position', [scrsz(3)/3 scrsz(4) scrsz(3)/3 scrsz(4)/2.5]);
        plot(tme, xcorr(PreMean - mean(PreMean), PostMean - mean(PostMean), MaxLag, 'coeff'), 'LineWidth', 2);
        axis tight
        title(SampleEpoch.protocolSettings.acquirinoCellBasename);
        LegendString = strcat('Pre = ', SampleEpoch.protocolSettings.PreSynapticHold);
        LegendString = strcat(LegendString, strcat('/Post = ', SampleEpoch.protocolSettings.PostSynapticHold));
        legend(LegendString, 'Location', 'Best');
        legend boxoff
        FileName = strcat(strcat(SampleEpoch.protocolSettings.acquirinoCellBasename, 'acor'), num2str(Leaves{cond}(leaf)));
        saveas(gcf, FileName, 'pdf');    

        figure(3); clf
        set(3, 'Position', [2*scrsz(3)/3 scrsz(4) scrsz(3)/3 scrsz(4)/2.5]);
        plot(tme, CrossCorrelation, 'LineWidth', 2);
        xlim([-0.1 0.1]);
        xlabel('time (sec)');
        ylabel('Correlation');
        title(SampleEpoch.protocolSettings.acquirinoCellBasename);
        legend(LegendString, 'Location', 'Best');
        legend boxoff
        FileName = strcat(strcat(SampleEpoch.protocolSettings.acquirinoCellBasename, 'xcor'), num2str(Leaves{cond}(leaf)));
        saveas(gcf, FileName, 'pdf');    
        pause(1);
    end 
end


%%