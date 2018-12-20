function[] = frequency_analysis_plots(rhods, thetads, angleds, datarun002, ds,StimComb)

maxVector = [0 1];
concatVector  =[0 1];
path = 0;
binSize = 0.1;
spPerd = unique(StimComb(:,1));
tempPerd = unique(StimComb(:,2));
maxStr(1,:) = 'Max';
maxStr(2,:) = 'Sum';

for a = 1:length(spPerd)
    for b = 1:length(tempPerd)
        [sigPowerPeakFreqAveMax, f1MagAveMax, f2MagAveMax, f2f1RatioAveMax] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 0, StimComb, spPerd(a,1), tempPerd(b,1), path, binSize);
        [sigPowerPeakFreqConMax, f1MagConMax, f2MagConMax, f2f1RatioConMax] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 1, StimComb, spPerd(a,1), tempPerd(b,1), path, binSize);
        [sigPowerPeakFreqAveSum, f1MagAveSum, f2MagAveSum, f2f1RatioAveSum] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 0, StimComb, spPerd(a,1), tempPerd(b,1), path, binSize);
        [sigPowerPeakFreqConSum, f1MagConSum, f2MagConSum, f2f1RatioConSum] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 1, StimComb, spPerd(a,1), tempPerd(b,1), path, binSize);
   
        subplot(4,4,1);
        plot(f1MagAveMax, f2MagAveMax, 'o');
        xlabel('F1 - AveSpMaxFR');
        ylabel('F2 - AveSpMaxFR');
        subplot(4,4,5);
        plot(f1MagConMax, f2MagConMax, 'o');
        xlabel('F1 - ConSpMaxFR');
        ylabel('F2 - ConSpMaxFR');
        subplot(4,4,9);
        plot(f1MagAveSum, f2MagAveSum, 'o');
        xlabel('F1 - AveSpSumFR');
        ylabel('F2 - AveSpSumFR');
        subplot(4,4,13);
        plot(f1MagConSum, f2MagConSum, 'o');
        xlabel('F1 - ConSpSumFR');
        ylabel('F2 - ConSpSumFR');
            
        subplot(4,4,2);
        hist(f2f1RatioAveMax, 40);
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,6);
        hist(f2f1RatioConMax, 40);
        xlabel('F2/F1 - concatSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,10);
        hist(f2f1RatioAveSum, 40);
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('No of observations');
        subplot(4,4,14);
        hist(f2f1RatioConSum, 40);
        xlabel('F2/F1 - concatSpSumFR');
        ylabel('No of observations');
            
        subplot(4,4,3);
        plot(f2f1RatioAveMax, f2f1RatioConMax, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - concatSpMaxFR');
        subplot(4,4,7);
        plot(f2f1RatioAveMax, f2f1RatioAveSum, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - aveSpSumFR');
        subplot(4,4,11);
        plot(f2f1RatioAveSum, f2f1RatioConSum, 'o');
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('F2/F1 - concatSpSumFR');
        subplot(4,4,15);
        plot(f2f1RatioConMax, f2f1RatioConSum, 'o');
        xlabel('F2/F1 - conSpMaxFR');
        ylabel('F2/F1 - conSpSumFR');
        
        subplot(4,4,4);
        plot(1:length(f2f1RatioConMax), f2f1RatioConMax, 'o');
        hold on;
        plot(1:length(f2f1RatioAveMax), f2f1RatioAveMax, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveMax red - F2/F1ConcatMax-blue');
        subplot(4,4,8)
        plot(1:length(f2f1RatioAveMax), f2f1RatioAveMax, 'o');
        hold on;
        plot(1:length(f2f1RatioAveSum), f2f1RatioAveSum, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveSum red - F2/F1AveMax-blue');
        subplot(4,4,12)
        plot(1:length(f2f1RatioAveSum), f2f1RatioAveSum, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum), f2f1RatioConSum, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1AveSum-blue');
        subplot(4,4,16)
        plot(1:length(f2f1RatioConMax), f2f1RatioConMax, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum), f2f1RatioConSum, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1ConMax-blue');
          
         pause;                

        end
    end
end

[sigPowerPeakFreqAveMax64, f1MagAveMax64, f2MagAveMax64, f2f1RatioAveMax64] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 0, StimComb, 64, 64, path, binSize);
[sigPowerPeakFreqConMax64, f1MagConMax64, f2MagConMax64, f2f1RatioConMax64] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 1, StimComb, 64, 64, path, binSize);
[sigPowerPeakFreqAveSum64, f1MagAveSum64, f2MagAveSum64, f2f1RatioAveSum64] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 0, StimComb, 64, 64, path, binSize);
[sigPowerPeakFreqConSum64, f1MagConSum64, f2MagConSum64, f2f1RatioConSum64] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 1, StimComb, 64, 64, path, binSize);

[sigPowerPeakFreqAveMax256, f1MagAveMax256, f2MagAveMax256, f2f1RatioAveMax256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 0, StimComb, 64, 256, path, binSize);
[sigPowerPeakFreqConMax256, f1MagConMax256, f2MagConMax256, f2f1RatioConMax256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 1, StimComb, 64, 256, path, binSize);
[sigPowerPeakFreqAveSum256, f1MagAveSum256, f2MagAveSum256, f2f1RatioAveSum256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 0, StimComb, 64, 256, path, binSize);
[sigPowerPeakFreqConSum256, f1MagConSum256, f2MagConSum256, f2f1RatioConSum256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 1, StimComb, 64, 256, path, binSize);

[sigPowerPeakFreqAveMax512, f1MagAveMax512, f2MagAveMax512, f2f1RatioAveMax512] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 0, StimComb, 64, 512, path, binSize);
[sigPowerPeakFreqConMax512, f1MagConMax512, f2MagConMax512, f2f1RatioConMax512] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 1, StimComb, 64, 512, path, binSize);
[sigPowerPeakFreqAveSum512, f1MagAveSum512, f2MagAveSum512, f2f1RatioAveSum512] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 0, StimComb, 64, 512, path, binSize);
[sigPowerPeakFreqConSum512, f1MagConSum512, f2MagConSum512, f2f1RatioConSum512] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 1, StimComb, 64, 512, path, binSize);


%%%%%% 64
        subplot(4,4,1);
        plot(f1MagAveMax64, f2MagAveMax64, 'o');
        xlabel('F1 - AveSpMaxFR');
        ylabel('F2 - AveSpMaxFR');
        subplot(4,4,5);
        plot(f1MagConMax64, f2MagConMax64, 'o');
        xlabel('F1 - ConSpMaxFR');
        ylabel('F2 - ConSpMaxFR');
        subplot(4,4,9);
        plot(f1MagAveSum64, f2MagAveSum64, 'o');
        xlabel('F1 - AveSpSumFR');
        ylabel('F2 - AveSpSumFR');
        subplot(4,4,13);
        plot(f1MagConSum64, f2MagConSum64, 'o');
        xlabel('F1 - ConSpSumFR');
        ylabel('F2 - ConSpSumFR');
            
        subplot(4,4,2);
        hist(f2f1RatioAveMax64, 40);
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,6);
        hist(f2f1RatioConMax64, 40);
        xlabel('F2/F1 - concatSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,10);
        hist(f2f1RatioAveSum64, 40);
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('No of observations');
        subplot(4,4,14);
        hist(f2f1RatioConSum64, 40);
        xlabel('F2/F1 - concatSpSumFR');
        ylabel('No of observations');
            
        subplot(4,4,3);
        plot(f2f1RatioAveMax64, f2f1RatioConMax64, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - concatSpMaxFR');
        subplot(4,4,7);
        plot(f2f1RatioAveMax64, f2f1RatioAveSum64, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - aveSpSumFR');
        subplot(4,4,11);
        plot(f2f1RatioAveSum64, f2f1RatioConSum64, 'o');
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('F2/F1 - concatSpSumFR');
        subplot(4,4,15);
        plot(f2f1RatioConMax64, f2f1RatioConSum64, 'o');
        xlabel('F2/F1 - conSpMaxFR');
        ylabel('F2/F1 - conSpSumFR');
        
        subplot(4,4,4);
        plot(1:length(f2f1RatioConMax64), f2f1RatioConMax64, 'o');
        hold on;
        plot(1:length(f2f1RatioAveMax64), f2f1RatioAveMax64, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveMax red - F2/F1ConcatMax-blue');
        subplot(4,4,8)
        plot(1:length(f2f1RatioAveMax64), f2f1RatioAveMax64, 'o');
        hold on;
        plot(1:length(f2f1RatioAveSum64), f2f1RatioAveSum64, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveSum red - F2/F1AveMax-blue');
        subplot(4,4,12)
        plot(1:length(f2f1RatioAveSum64), f2f1RatioAveSum64, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum64), f2f1RatioConSum64, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1AveSum-blue');
        subplot(4,4,16)
        plot(1:length(f2f1RatioConMax64), f2f1RatioConMax64, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum64), f2f1RatioConSum64, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1ConMax-blue');

%%%%%%%% 256

        subplot(4,4,1);
        plot(f1MagAveMax256, f2MagAveMax256, 'o');
        xlabel('F1 - AveSpMaxFR');
        ylabel('F2 - AveSpMaxFR');
        subplot(4,4,5);
        plot(f1MagConMax256, f2MagConMax256, 'o');
        xlabel('F1 - ConSpMaxFR');
        ylabel('F2 - ConSpMaxFR');
        subplot(4,4,9);
        plot(f1MagAveSum256, f2MagAveSum256, 'o');
        xlabel('F1 - AveSpSumFR');
        ylabel('F2 - AveSpSumFR');
        subplot(4,4,13);
        plot(f1MagConSum256, f2MagConSum256, 'o');
        xlabel('F1 - ConSpSumFR');
        ylabel('F2 - ConSpSumFR');
            
        subplot(4,4,2);
        hist(f2f1RatioAveMax256, 40);
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,6);
        hist(f2f1RatioConMax256, 40);
        xlabel('F2/F1 - concatSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,10);
        hist(f2f1RatioAveSum256, 40);
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('No of observations');
        subplot(4,4,14);
        hist(f2f1RatioConSum256, 40);
        xlabel('F2/F1 - concatSpSumFR');
        ylabel('No of observations');
            
        subplot(4,4,3);
        plot(f2f1RatioAveMax256, f2f1RatioConMax256, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - concatSpMaxFR');
        subplot(4,4,7);
        plot(f2f1RatioAveMax256, f2f1RatioAveSum256, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - aveSpSumFR');
        subplot(4,4,11);
        plot(f2f1RatioAveSum256, f2f1RatioConSum256, 'o');
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('F2/F1 - concatSpSumFR');
        subplot(4,4,15);
        plot(f2f1RatioConMax256, f2f1RatioConSum256, 'o');
        xlabel('F2/F1 - conSpMaxFR');
        ylabel('F2/F1 - conSpSumFR');
        
        subplot(4,4,4);
        plot(1:length(f2f1RatioConMax256), f2f1RatioConMax256, 'o');
        hold on;
        plot(1:length(f2f1RatioAveMax256), f2f1RatioAveMax256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveMax red - F2/F1ConcatMax-blue');
        subplot(4,4,8)
        plot(1:length(f2f1RatioAveMax256), f2f1RatioAveMax256, 'o');
        hold on;
        plot(1:length(f2f1RatioAveSum256), f2f1RatioAveSum256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveSum red - F2/F1AveMax-blue');
        subplot(4,4,12)
        plot(1:length(f2f1RatioAveSum256), f2f1RatioAveSum256, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum256), f2f1RatioConSum256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1AveSum-blue');
        subplot(4,4,16)
        plot(1:length(f2f1RatioConMax256), f2f1RatioConMax256, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum256), f2f1RatioConSum256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1ConMax-blue');
        
        
%%%%%%%%%%% 512
        subplot(4,4,1);
        plot(f1MagAveMax512, f2MagAveMax512, 'o');
        xlabel('F1 - AveSpMaxFR');
        ylabel('F2 - AveSpMaxFR');
        subplot(4,4,5);
        plot(f1MagConMax512, f2MagConMax512, 'o');
        xlabel('F1 - ConSpMaxFR');
        ylabel('F2 - ConSpMaxFR');
        subplot(4,4,9);
        plot(f1MagAveSum512, f2MagAveSum512, 'o');
        xlabel('F1 - AveSpSumFR');
        ylabel('F2 - AveSpSumFR');
        subplot(4,4,13);
        plot(f1MagConSum512, f2MagConSum512, 'o');
        xlabel('F1 - ConSpSumFR');
        ylabel('F2 - ConSpSumFR');
            
        subplot(4,4,2);
        hist(f2f1RatioAveMax512, 40);
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,6);
        hist(f2f1RatioConMax512, 40);
        xlabel('F2/F1 - concatSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,10);
        hist(f2f1RatioAveSum512, 40);
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('No of observations');
        subplot(4,4,14);
        hist(f2f1RatioConSum512, 40);
        xlabel('F2/F1 - concatSpSumFR');
        ylabel('No of observations');
            
        subplot(4,4,3);
        plot(f2f1RatioAveMax512, f2f1RatioConMax512, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - concatSpMaxFR');
        subplot(4,4,7);
        plot(f2f1RatioAveMax512, f2f1RatioAveSum512, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - aveSpSumFR');
        subplot(4,4,11);
        plot(f2f1RatioAveSum512, f2f1RatioConSum512, 'o');
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('F2/F1 - concatSpSumFR');
        subplot(4,4,15);
        plot(f2f1RatioConMax512, f2f1RatioConSum512, 'o');
        xlabel('F2/F1 - conSpMaxFR');
        ylabel('F2/F1 - conSpSumFR');
        
        subplot(4,4,4);
        plot(1:length(f2f1RatioConMax512), f2f1RatioConMax512, 'o');
        hold on;
        plot(1:length(f2f1RatioAveMax512), f2f1RatioAveMax512, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveMax red - F2/F1ConcatMax-blue');
        subplot(4,4,8)
        plot(1:length(f2f1RatioAveMax512), f2f1RatioAveMax512, 'o');
        hold on;
        plot(1:length(f2f1RatioAveSum512), f2f1RatioAveSum512, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveSum red - F2/F1AveMax-blue');
        subplot(4,4,12)
        plot(1:length(f2f1RatioAveSum512), f2f1RatioAveSum512, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum512), f2f1RatioConSum512, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1AveSum-blue');
        subplot(4,4,16)
        plot(1:length(f2f1RatioConMax512), f2f1RatioConMax512, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum512), f2f1RatioConSum512, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1ConMax-blue');

        
        
        %%%%%%%%%speeds
        plot(f2f1RatioAveMax64, f2f1RatioAveMax256, 'o');
        


%%%%%%%%%%

[datarun000] = load_dsdata('/Analysis/sravi/2012-10-31-0/data000-1800-7200/', 'data000-map/data000-map', 0, 0, 1);
[datarun001] = load_dsdata('/Analysis/sravi/2012-10-31-0/data000-1800-7200/', 'data001-map/data001-map', 0, 0, 0);
[datarun002] = load_dsdata('/Analysis/sravi/2012-10-31-0/data000-1800-7200/', 'data002-map/data002-map', 1, '/stimuli/s02', 0);
cellids = intersect(intersect(datarun000.cell_ids, datarun001.cell_ids),datarun002.cell_ids);
[NumSpikesCell, StimComb] = get_spikescellstim(datarun002, datarun002.cell_ids, 0);
NumSpikesCellIn = NumSpikesCell(get_cell_indices(datarun002, cellids), :);
[magin dsindexin magmaxin magavein anglein] = dscellanalysis(NumSpikesCellIn, StimComb);
plot(magin{1,1}, magin{2,1},'o');
plot(magin{1,2}, magin{2,2},'o');

ginput(1);
Y = get(get(gca,'Children'),'YData');
X = get(get(gca,'Children'),'XData');
h = impoly;
accepted_pos = wait(h);
hold on;
drawPolygon(accepted_pos(:,1), accepted_pos(:,2));
in = inpolygon(X, Y, accepted_pos(:,1), accepted_pos(:,2));
plot(X(in),Y(in),'r+');
cid = 1:length(X);
chos = cid(in);
ds = cellids(1, chos);
NumSpikesCellDS = NumSpikesCell(get_cell_indices(datarun002, ds), :);
[magds dsindexds magmaxds magaveds angleds rhods thetads numds] = dscellanalysis(NumSpikesCellDS, StimComb);



%%%%%%%%%%%%%%

[datarun000] = load_dsdata('/Analysis/sravi/2012-10-10-1/data000-1800-7200/', 'data000-map/data000-map', 0, 0, 1);
[datarun001] = load_dsdata('/Analysis/sravi/2012-10-10-1/data000-1800-7200/', 'data001-map/data001-map', 0, 0, 0);
[datarun002] = load_dsdata('/Analysis/sravi/2012-10-10-1/data000-1800-7200/', 'data002-map/data002-map', 1, '/stimuli/s02', 0);

%%%%%%

[sigPowerPeakFreqAveMax64, f1MagAveMax64, f2MagAveMax64, f2f1RatioAveMax64] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 0, StimComb, 256, 64, path, binSize);
[sigPowerPeakFreqConMax64, f1MagConMax64, f2MagConMax64, f2f1RatioConMax64] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 1, StimComb, 256, 64, path, binSize);
[sigPowerPeakFreqAveSum64, f1MagAveSum64, f2MagAveSum64, f2f1RatioAveSum64] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 0, StimComb, 256, 64, path, binSize);
[sigPowerPeakFreqConSum64, f1MagConSum64, f2MagConSum64, f2f1RatioConSum64] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 1, StimComb, 256, 64, path, binSize);

[sigPowerPeakFreqAveMax256, f1MagAveMax256, f2MagAveMax256, f2f1RatioAveMax256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 0, StimComb, 256, 256, path, binSize);
[sigPowerPeakFreqConMax256, f1MagConMax256, f2MagConMax256, f2f1RatioConMax256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 1, StimComb, 256, 256, path, binSize);
[sigPowerPeakFreqAveSum256, f1MagAveSum256, f2MagAveSum256, f2f1RatioAveSum256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 0, StimComb, 256, 256, path, binSize);
[sigPowerPeakFreqConSum256, f1MagConSum256, f2MagConSum256, f2f1RatioConSum256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 1, StimComb, 256, 256, path, binSize);


 subplot(4,4,1);
        plot(f1MagAveMax64, f2MagAveMax64, 'o');
        xlabel('F1 - AveSpMaxFR');
        ylabel('F2 - AveSpMaxFR');
        subplot(4,4,5);
        plot(f1MagConMax64, f2MagConMax64, 'o');
        xlabel('F1 - ConSpMaxFR');
        ylabel('F2 - ConSpMaxFR');
        subplot(4,4,9);
        plot(f1MagAveSum64, f2MagAveSum64, 'o');
        xlabel('F1 - AveSpSumFR');
        ylabel('F2 - AveSpSumFR');
        subplot(4,4,13);
        plot(f1MagConSum64, f2MagConSum64, 'o');
        xlabel('F1 - ConSpSumFR');
        ylabel('F2 - ConSpSumFR');
            
        subplot(4,4,2);
        hist(f2f1RatioAveMax64, 40);
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,6);
        hist(f2f1RatioConMax64, 40);
        xlabel('F2/F1 - concatSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,10);
        hist(f2f1RatioAveSum64, 40);
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('No of observations');
        subplot(4,4,14);
        hist(f2f1RatioConSum64, 40);
        xlabel('F2/F1 - concatSpSumFR');
        ylabel('No of observations');
            
        subplot(4,4,3);
        plot(f2f1RatioAveMax64, f2f1RatioConMax64, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - concatSpMaxFR');
        subplot(4,4,7);
        plot(f2f1RatioAveMax64, f2f1RatioAveSum64, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - aveSpSumFR');
        subplot(4,4,11);
        plot(f2f1RatioAveSum64, f2f1RatioConSum64, 'o');
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('F2/F1 - concatSpSumFR');
        subplot(4,4,15);
        plot(f2f1RatioConMax64, f2f1RatioConSum64, 'o');
        xlabel('F2/F1 - conSpMaxFR');
        ylabel('F2/F1 - conSpSumFR');
        
        subplot(4,4,4);
        plot(1:length(f2f1RatioConMax64), f2f1RatioConMax64, 'o');
        hold on;
        plot(1:length(f2f1RatioAveMax64), f2f1RatioAveMax64, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveMax red - F2/F1ConcatMax-blue');
        subplot(4,4,8)
        plot(1:length(f2f1RatioAveMax64), f2f1RatioAveMax64, 'o');
        hold on;
        plot(1:length(f2f1RatioAveSum64), f2f1RatioAveSum64, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveSum red - F2/F1AveMax-blue');
        subplot(4,4,12)
        plot(1:length(f2f1RatioAveSum64), f2f1RatioAveSum64, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum64), f2f1RatioConSum64, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1AveSum-blue');
        subplot(4,4,16)
        plot(1:length(f2f1RatioConMax64), f2f1RatioConMax64, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum64), f2f1RatioConSum64, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1ConMax-blue');

%%%%%%%% 256

        subplot(4,4,1);
        plot(f1MagAveMax256, f2MagAveMax256, 'o');
        xlabel('F1 - AveSpMaxFR');
        ylabel('F2 - AveSpMaxFR');
        subplot(4,4,5);
        plot(f1MagConMax256, f2MagConMax256, 'o');
        xlabel('F1 - ConSpMaxFR');
        ylabel('F2 - ConSpMaxFR');
        subplot(4,4,9);
        plot(f1MagAveSum256, f2MagAveSum256, 'o');
        xlabel('F1 - AveSpSumFR');
        ylabel('F2 - AveSpSumFR');
        subplot(4,4,13);
        plot(f1MagConSum256, f2MagConSum256, 'o');
        xlabel('F1 - ConSpSumFR');
        ylabel('F2 - ConSpSumFR');
            
        subplot(4,4,2);
        hist(f2f1RatioAveMax256, 40);
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,6);
        hist(f2f1RatioConMax256, 40);
        xlabel('F2/F1 - concatSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,10);
        hist(f2f1RatioAveSum256, 40);
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('No of observations');
        subplot(4,4,14);
        hist(f2f1RatioConSum256, 40);
        xlabel('F2/F1 - concatSpSumFR');
        ylabel('No of observations');
            
        subplot(4,4,3);
        plot(f2f1RatioAveMax256, f2f1RatioConMax256, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - concatSpMaxFR');
        subplot(4,4,7);
        plot(f2f1RatioAveMax256, f2f1RatioAveSum256, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - aveSpSumFR');
        subplot(4,4,11);
        plot(f2f1RatioAveSum256, f2f1RatioConSum256, 'o');
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('F2/F1 - concatSpSumFR');
        subplot(4,4,15);
        plot(f2f1RatioConMax256, f2f1RatioConSum256, 'o');
        xlabel('F2/F1 - conSpMaxFR');
        ylabel('F2/F1 - conSpSumFR');
        
        subplot(4,4,4);
        plot(1:length(f2f1RatioConMax256), f2f1RatioConMax256, 'o');
        hold on;
        plot(1:length(f2f1RatioAveMax256), f2f1RatioAveMax256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveMax red - F2/F1ConcatMax-blue');
        subplot(4,4,8)
        plot(1:length(f2f1RatioAveMax256), f2f1RatioAveMax256, 'o');
        hold on;
        plot(1:length(f2f1RatioAveSum256), f2f1RatioAveSum256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveSum red - F2/F1AveMax-blue');
        subplot(4,4,12)
        plot(1:length(f2f1RatioAveSum256), f2f1RatioAveSum256, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum256), f2f1RatioConSum256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1AveSum-blue');
        subplot(4,4,16)
        plot(1:length(f2f1RatioConMax256), f2f1RatioConMax256, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum256), f2f1RatioConSum256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1ConMax-blue');
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%2012-10-31-0%%%%%%%


figure(5)
plot(magin{1,1}, magin{2,1},'o');
plot(magin{1,2}, magin{2,2},'o');
plot(magin{1,1}, magin{1,2},'o');
plot(magin{2,1}, magin{2,2},'o');

ginput(1);
Y = get(get(gca,'Children'),'YData');
X = get(get(gca,'Children'),'XData');
h = impoly;
accepted_pos = wait(h);
hold on;
drawPolygon(accepted_pos(:,1), accepted_pos(:,2));
in = inpolygon(X, Y, accepted_pos(:,1), accepted_pos(:,2));
plot(X(in),Y(in),'r+');
cid = 1:length(X);
chos = cid(in);
ds = cellids(1, chos);

NumSpikesCellDS = NumSpikesCell(get_cell_indices(datarun002, ds), :);
[magds dsindexds magmaxds magaveds angleds rhods thetads numds] = dscellanalysis(NumSpikesCellDS, StimComb);

maxVal = 1;
concat = 1;
path = 0;
binSize = 0.1;

[sigPowerPeakFreq6464, f1Mag6464, f2Mag6464, f2f1Ratio6464] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, maxVal, concat, StimComb, 64, 64, path, binSize);
[sigPowerPeakFreq64256, f1Mag64256, f2Mag64256, f2f1Ratio64256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, maxVal, concat, StimComb, 64, 256, path, binSize);

[sigPowerPeakFreq25664, f1Mag25664, f2Mag25664, f2f1Ratio25664] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, maxVal, concat, StimComb, 256, 64, path, binSize);
[sigPowerPeakFreq256256, f1Mag256256, f2Mag256256, f2f1Ratio256256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, maxVal, concat, StimComb, 256, 256, path, binSize);



[sigPowerPeakFreq, f1Mag, f2Mag, f2f1Ratio] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 1, StimComb, 512, 128, 0, 0.1);


  
subplot(4, 2, 1)
plot(f1Mag6464, f2Mag6464, 'o');
xlabel('F1 - ConSpMaxFR');
ylabel('F2 - ConSpMaxFR');
title('64 - 64');
   
subplot(4, 2, 2)
hist(f2f1Ratio6464, 40);
xlabel('F2/F1 - ConSpMaxFR');
ylabel('No of observations');
title('64 - 64');
       
subplot(4, 2, 3)
plot(f1Mag64256, f2Mag64256, 'o');
xlabel('F1 - ConSpMaxFR');
ylabel('F2 - ConSpMaxFR');
title('64 - 256');
   
subplot(4, 2, 4)
hist(f2f1Ratio64256, 40);
xlabel('F2/F1 - ConSpMaxFR');
ylabel('No of observations');
title('64 - 256');
       
subplot(4, 2, 5)
plot(f1Mag25664, f2Mag25664, 'o');
xlabel('F1 - ConSpMaxFR');
ylabel('F2 - ConSpMaxFR');
title('256 - 64');
   
subplot(4, 2, 6)
hist(f2f1Ratio25664, 40);
xlabel('F2/F1 - ConSpMaxFR');
ylabel('No of observations');
title('256 - 64');
       
subplot(4, 2, 7)
plot(f1Mag256256, f2Mag256256, 'o');
xlabel('F1 - ConSpMaxFR');
ylabel('F2 - ConSpMaxFR');
title('256 - 256');
   
subplot(4, 2, 8)
hist(f2f1Ratio256256, 40);
xlabel('F2/F1 - ConSpMaxFR');
ylabel('No of observations');
title('256 - 256');

%same spatial diff temporal
plot(f2f1Ratio6464, f2f1Ratio64256, 'o');
plot(f2f1Ratio25664, f2f1Ratio256256, 'o');
%same temp diff spatial
plot(f2f1Ratio6464, f2f1Ratio25664, 'o');
plot(f2f1Ratio64256, f2f1Ratio256256, 'o');



%%%%%

[sigPowerPeakFreqAveMax6464, f1MagAveMax6464, f2MagAveMax6464, f2f1RatioAveMax6464] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 0, StimComb, 64, 64, path, binSize);
[sigPowerPeakFreqConMax6464, f1MagConMax6464, f2MagConMax6464, f2f1RatioConMax6464] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 1, StimComb, 64, 64, path, binSize);
[sigPowerPeakFreqAveSum6464, f1MagAveSum6464, f2MagAveSum6464, f2f1RatioAveSum6464] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 0, StimComb, 64, 64, path, binSize);
[sigPowerPeakFreqConSum6464, f1MagConSum6464, f2MagConSum6464, f2f1RatioConSum6464] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 1, StimComb, 64, 64, path, binSize);

[sigPowerPeakFreqAveMax64256, f1MagAveMax64256, f2MagAveMax64256, f2f1RatioAveMax64256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 0, StimComb, 64, 256, path, binSize);
[sigPowerPeakFreqConMax64256, f1MagConMax64256, f2MagConMax64256, f2f1RatioConMax64256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 1, StimComb, 64, 256, path, binSize);
[sigPowerPeakFreqAveSum64256, f1MagAveSum64256, f2MagAveSum64256, f2f1RatioAveSum64256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 0, StimComb, 64, 256, path, binSize);
[sigPowerPeakFreqConSum64256, f1MagConSum64256, f2MagConSum64256, f2f1RatioConSum64256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 1, StimComb, 64, 256, path, binSize);

[sigPowerPeakFreqAveMax25664, f1MagAveMax25664, f2MagAveMax25664, f2f1RatioAveMax25664] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 0, StimComb, 256, 64, path, binSize);
[sigPowerPeakFreqConMax25664, f1MagConMax25664, f2MagConMax25664, f2f1RatioConMax25664] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 1, StimComb, 256, 64, path, binSize);
[sigPowerPeakFreqAveSum25664, f1MagAveSum25664, f2MagAveSum25664, f2f1RatioAveSum25664] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 0, StimComb, 256, 64, path, binSize);
[sigPowerPeakFreqConSum25664, f1MagConSum25664, f2MagConSum25664, f2f1RatioConSum25664] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 1, StimComb, 256, 64, path, binSize);

[sigPowerPeakFreqAveMax256256, f1MagAveMax256256, f2MagAveMax256256, f2f1RatioAveMax256256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 0, StimComb, 256, 256, path, binSize);
[sigPowerPeakFreqConMax256256, f1MagConMax256256, f2MagConMax256256, f2f1RatioConMax256256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, 1, StimComb, 256, 256, path, binSize);
[sigPowerPeakFreqAveSum256256, f1MagAveSum256256, f2MagAveSum256256, f2f1RatioAveSum256256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 0, StimComb, 256, 256, path, binSize);
[sigPowerPeakFreqConSum256256, f1MagConSum256256, f2MagConSum256256, f2f1RatioConSum256256] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, 1, StimComb, 256, 256, path, binSize);



figure(1);
 subplot(4,4,1);
        plot(f1MagAveMax6464, f2MagAveMax6464, 'o');
        xlabel('F1 - AveSpMaxFR');
        ylabel('F2 - AveSpMaxFR');
        subplot(4,4,5);
        plot(f1MagConMax6464, f2MagConMax6464, 'o');
        xlabel('F1 - ConSpMaxFR');
        ylabel('F2 - ConSpMaxFR');
        subplot(4,4,9);
        plot(f1MagAveSum6464, f2MagAveSum6464, 'o');
        xlabel('F1 - AveSpSumFR');
        ylabel('F2 - AveSpSumFR');
        subplot(4,4,13);
        plot(f1MagConSum6464, f2MagConSum6464, 'o');
        xlabel('F1 - ConSpSumFR');
        ylabel('F2 - ConSpSumFR');
            
        subplot(4,4,2);
        hist(f2f1RatioAveMax6464, 40);
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,6);
        hist(f2f1RatioConMax6464, 40);
        xlabel('F2/F1 - concatSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,10);
        hist(f2f1RatioAveSum6464, 40);
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('No of observations');
        subplot(4,4,14);
        hist(f2f1RatioConSum6464, 40);
        xlabel('F2/F1 - concatSpSumFR');
        ylabel('No of observations');
            
        subplot(4,4,3);
        plot(f2f1RatioAveMax6464, f2f1RatioConMax6464, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - concatSpMaxFR');
        subplot(4,4,7);
        plot(f2f1RatioAveMax6464, f2f1RatioAveSum6464, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - aveSpSumFR');
        subplot(4,4,11);
        plot(f2f1RatioAveSum6464, f2f1RatioConSum6464, 'o');
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('F2/F1 - concatSpSumFR');
        subplot(4,4,15);
        plot(f2f1RatioConMax6464, f2f1RatioConSum6464, 'o');
        xlabel('F2/F1 - conSpMaxFR');
        ylabel('F2/F1 - conSpSumFR');
        
        subplot(4,4,4);
        plot(1:length(f2f1RatioConMax6464), f2f1RatioConMax6464, 'o');
        hold on;
        plot(1:length(f2f1RatioAveMax6464), f2f1RatioAveMax6464, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveMax red - F2/F1ConcatMax-blue');
        subplot(4,4,8)
        plot(1:length(f2f1RatioAveMax6464), f2f1RatioAveMax6464, 'o');
        hold on;
        plot(1:length(f2f1RatioAveSum6464), f2f1RatioAveSum6464, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveSum red - F2/F1AveMax-blue');
        subplot(4,4,12)
        plot(1:length(f2f1RatioAveSum6464), f2f1RatioAveSum6464, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum6464), f2f1RatioConSum6464, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1AveSum-blue');
        subplot(4,4,16)
        plot(1:length(f2f1RatioConMax6464), f2f1RatioConMax6464, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum6464), f2f1RatioConSum6464, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1ConMax-blue');

%%%%%%%% 256
figure(2);
        subplot(4,4,1);
        plot(f1MagAveMax64256, f2MagAveMax64256, 'o');
        xlabel('F1 - AveSpMaxFR');
        ylabel('F2 - AveSpMaxFR');
        subplot(4,4,5);
        plot(f1MagConMax64256, f2MagConMax64256, 'o');
        xlabel('F1 - ConSpMaxFR');
        ylabel('F2 - ConSpMaxFR');
        subplot(4,4,9);
        plot(f1MagAveSum64256, f2MagAveSum64256, 'o');
        xlabel('F1 - AveSpSumFR');
        ylabel('F2 - AveSpSumFR');
        subplot(4,4,13);
        plot(f1MagConSum64256, f2MagConSum64256, 'o');
        xlabel('F1 - ConSpSumFR');
        ylabel('F2 - ConSpSumFR');
            
        subplot(4,4,2);
        hist(f2f1RatioAveMax64256, 40);
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,6);
        hist(f2f1RatioConMax64256, 40);
        xlabel('F2/F1 - concatSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,10);
        hist(f2f1RatioAveSum64256, 40);
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('No of observations');
        subplot(4,4,14);
        hist(f2f1RatioConSum64256, 40);
        xlabel('F2/F1 - concatSpSumFR');
        ylabel('No of observations');
            
        subplot(4,4,3);
        plot(f2f1RatioAveMax64256, f2f1RatioConMax64256, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - concatSpMaxFR');
        subplot(4,4,7);
        plot(f2f1RatioAveMax64256, f2f1RatioAveSum64256, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - aveSpSumFR');
        subplot(4,4,11);
        plot(f2f1RatioAveSum64256, f2f1RatioConSum64256, 'o');
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('F2/F1 - concatSpSumFR');
        subplot(4,4,15);
        plot(f2f1RatioConMax64256, f2f1RatioConSum64256, 'o');
        xlabel('F2/F1 - conSpMaxFR');
        ylabel('F2/F1 - conSpSumFR');
        
        subplot(4,4,4);
        plot(1:length(f2f1RatioConMax64256), f2f1RatioConMax64256, 'o');
        hold on;
        plot(1:length(f2f1RatioAveMax64256), f2f1RatioAveMax64256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveMax red - F2/F1ConcatMax-blue');
        subplot(4,4,8)
        plot(1:length(f2f1RatioAveMax64256), f2f1RatioAveMax64256, 'o');
        hold on;
        plot(1:length(f2f1RatioAveSum64256), f2f1RatioAveSum64256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveSum red - F2/F1AveMax-blue');
        subplot(4,4,12)
        plot(1:length(f2f1RatioAveSum64256), f2f1RatioAveSum64256, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum64256), f2f1RatioConSum64256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1AveSum-blue');
        subplot(4,4,16)
        plot(1:length(f2f1RatioConMax64256), f2f1RatioConMax64256, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum64256), f2f1RatioConSum64256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1ConMax-blue');
        



%%%%spatial = 256 temp = both


figure(3);
subplot(4,4,1);
        plot(f1MagAveMax25664, f2MagAveMax25664, 'o');
        xlabel('F1 - AveSpMaxFR');
        ylabel('F2 - AveSpMaxFR');
        subplot(4,4,5);
        plot(f1MagConMax25664, f2MagConMax25664, 'o');
        xlabel('F1 - ConSpMaxFR');
        ylabel('F2 - ConSpMaxFR');
        subplot(4,4,9);
        plot(f1MagAveSum25664, f2MagAveSum25664, 'o');
        xlabel('F1 - AveSpSumFR');
        ylabel('F2 - AveSpSumFR');
        subplot(4,4,13);
        plot(f1MagConSum25664, f2MagConSum25664, 'o');
        xlabel('F1 - ConSpSumFR');
        ylabel('F2 - ConSpSumFR');
            
        subplot(4,4,2);
        hist(f2f1RatioAveMax25664, 40);
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,6);
        hist(f2f1RatioConMax25664, 40);
        xlabel('F2/F1 - concatSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,10);
        hist(f2f1RatioAveSum25664, 40);
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('No of observations');
        subplot(4,4,14);
        hist(f2f1RatioConSum25664, 40);
        xlabel('F2/F1 - concatSpSumFR');
        ylabel('No of observations');
            
        subplot(4,4,3);
        plot(f2f1RatioAveMax25664, f2f1RatioConMax25664, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - concatSpMaxFR');
        subplot(4,4,7);
        plot(f2f1RatioAveMax25664, f2f1RatioAveSum25664, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - aveSpSumFR');
        subplot(4,4,11);
        plot(f2f1RatioAveSum25664, f2f1RatioConSum25664, 'o');
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('F2/F1 - concatSpSumFR');
        subplot(4,4,15);
        plot(f2f1RatioConMax25664, f2f1RatioConSum25664, 'o');
        xlabel('F2/F1 - conSpMaxFR');
        ylabel('F2/F1 - conSpSumFR');
        
        subplot(4,4,4);
        plot(1:length(f2f1RatioConMax25664), f2f1RatioConMax25664, 'o');
        hold on;
        plot(1:length(f2f1RatioAveMax25664), f2f1RatioAveMax25664, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveMax red - F2/F1ConcatMax-blue');
        subplot(4,4,8)
        plot(1:length(f2f1RatioAveMax25664), f2f1RatioAveMax25664, 'o');
        hold on;
        plot(1:length(f2f1RatioAveSum25664), f2f1RatioAveSum25664, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveSum red - F2/F1AveMax-blue');
        subplot(4,4,12)
        plot(1:length(f2f1RatioAveSum25664), f2f1RatioAveSum25664, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum25664), f2f1RatioConSum25664, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1AveSum-blue');
        subplot(4,4,16)
        plot(1:length(f2f1RatioConMax25664), f2f1RatioConMax25664, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum25664), f2f1RatioConSum25664, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1ConMax-blue');

%%%%%%%% 256
figure(4);
        subplot(4,4,1);
        plot(f1MagAveMax256256, f2MagAveMax256256, 'o');
        xlabel('F1 - AveSpMaxFR');
        ylabel('F2 - AveSpMaxFR');
        subplot(4,4,5);
        plot(f1MagConMax256256, f2MagConMax256256, 'o');
        xlabel('F1 - ConSpMaxFR');
        ylabel('F2 - ConSpMaxFR');
        subplot(4,4,9);
        plot(f1MagAveSum256256, f2MagAveSum256256, 'o');
        xlabel('F1 - AveSpSumFR');
        ylabel('F2 - AveSpSumFR');
        subplot(4,4,13);
        plot(f1MagConSum256256, f2MagConSum256256, 'o');
        xlabel('F1 - ConSpSumFR');
        ylabel('F2 - ConSpSumFR');
            
        subplot(4,4,2);
        hist(f2f1RatioAveMax256256, 40);
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,6);
        hist(f2f1RatioConMax256256, 40);
        xlabel('F2/F1 - concatSpMaxFR');
        ylabel('No of observations');
        subplot(4,4,10);
        hist(f2f1RatioAveSum256256, 40);
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('No of observations');
        subplot(4,4,14);
        hist(f2f1RatioConSum256256, 40);
        xlabel('F2/F1 - concatSpSumFR');
        ylabel('No of observations');
            
        subplot(4,4,3);
        plot(f2f1RatioAveMax256256, f2f1RatioConMax256256, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - concatSpMaxFR');
        subplot(4,4,7);
        plot(f2f1RatioAveMax256256, f2f1RatioAveSum256256, 'o');
        xlabel('F2/F1 - averageSpMaxFR');
        ylabel('F2/F1 - aveSpSumFR');
        subplot(4,4,11);
        plot(f2f1RatioAveSum256256, f2f1RatioConSum256256, 'o');
        xlabel('F2/F1 - averageSpSumFR');
        ylabel('F2/F1 - concatSpSumFR');
        subplot(4,4,15);
        plot(f2f1RatioConMax256256, f2f1RatioConSum256256, 'o');
        xlabel('F2/F1 - conSpMaxFR');
        ylabel('F2/F1 - conSpSumFR');
        
        subplot(4,4,4);
        plot(1:length(f2f1RatioConMax256256), f2f1RatioConMax256256, 'o');
        hold on;
        plot(1:length(f2f1RatioAveMax256256), f2f1RatioAveMax256256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveMax red - F2/F1ConcatMax-blue');
        subplot(4,4,8)
        plot(1:length(f2f1RatioAveMax256256), f2f1RatioAveMax256256, 'o');
        hold on;
        plot(1:length(f2f1RatioAveSum256256), f2f1RatioAveSum256256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1AveSum red - F2/F1AveMax-blue');
        subplot(4,4,12)
        plot(1:length(f2f1RatioAveSum256256), f2f1RatioAveSum256256, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum256256), f2f1RatioConSum256256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1AveSum-blue');
        subplot(4,4,16)
        plot(1:length(f2f1RatioConMax256256), f2f1RatioConMax256256, 'o');
        hold on;
        plot(1:length(f2f1RatioConSum256256), f2f1RatioConSum256256, 'o', 'Color','red');
        hold off;
        xlabel('cell');
        ylabel('F2/F1ConSum red - F2/F1ConMax-blue');
        






































%             titlechar = ['Sp ' num2str(spPerd(a,1)) 'Temp ' num2str(tempPerd(b,1)) 'V ' maxStr(c, :)];
%             title(titlechar);

    
% end   
    
%     sigPowerPeakFreq512max, f1Mag512max, f2Mag512max, f2f1Ratio512max] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 1, StimComb, 64, 512, path, binSize);
% 
% [sigPowerPeakFreq64mean, f1Mag64mean, f2Mag64mean, f2f1Ratio64mean] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, StimComb, 64, 64, path, binSize);
% [sigPowerPeakFreq256mean, f1Mag256mean, f2Mag256mean, f2f1Ratio256mean] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, StimComb, 64, 256, path, binSize);
% [sigPowerPeakFreq512mean, f1Mag512mean, f2Mag512mean, f2f1Ratio512mean] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, 0, StimComb, 64, 512, path, binSize);
