saveFigFlag = true ; % save plots to jpeg folder
saveFigPath = '/Users/jcafaro/Desktop/classificationFigs/' ;

% load data
NeuronStructPath = '/NeuronStruct' ;
load(NeuronStructPath)

numCells = length(NeuronStruct) ;
sampleRate = 20000 ;

% figures
CellFig = figure ; % figure of each cell

for cl = 1:numCells ;
    figure(CellFig); clf
    
    subplot(3,6,1:2)
    imagesc(NeuronStruct(cl).smallSquare.StvMax)
    hold on
    plot(NeuronStruct(cl).smallSquare.LikelyContour(1,:),NeuronStruct(cl).smallSquare.LikelyContour(2,:),'k')
    axis([0 400 0 400])
    title('small squares stv')
    
    subplot(3,6,3:4)
    imagesc(NeuronStruct(cl).bigSquare.StvMax)
    hold on
    plot(NeuronStruct(cl).bigSquare.Contours{1}(1,:),NeuronStruct(cl).bigSquare.Contours{1}(2,:),'k')
    axis([0 400 0 400])
    title('big squares stv')

    subplot(3,6,7:8)
    rasterPlot(NeuronStruct(cl).smallSquare.Spikes,[0,2])
    xlim([0,2])
    title('small squares trial with max response')
    
    subplot(3,6,9:10)
    rasterPlot(NeuronStruct(cl).bigSquare.Spikes,[0,2])
    xlim([0,2])
    title('big squares trial with max response')
    
    subplot(3,6,13:14)
    plot([1:length(NeuronStruct(cl).smallSquare.psth)]*1/sampleRate,NeuronStruct(cl).smallSquare.psth)
    xlim([0,2])
    ylim([0, max([NeuronStruct(cl).smallSquare.psth,NeuronStruct(cl).bigSquare.psth,NeuronStruct(cl).ffPulse.psth])])
    
    subplot(3,6,15:16)
    plot([1:length(NeuronStruct(cl).smallSquare.psth)]*1/sampleRate,NeuronStruct(cl).bigSquare.psth)
    xlim([0,2])
    ylim([0, max([NeuronStruct(cl).smallSquare.psth,NeuronStruct(cl).bigSquare.psth,NeuronStruct(cl).ffPulse.psth])])
    
    subplot(3,6,11:12)
    rasterPlot(NeuronStruct(cl).ffPulse.spikes,[0,2])
    xlim([0,2])
    title('full field step response')
    
    subplot(3,6,17:18)
    plot([1:length(NeuronStruct(cl).ffPulse.psth)]*1/sampleRate,NeuronStruct(cl).ffPulse.psth)
    xlim([0,2])
    ylim([0, max([NeuronStruct(cl).smallSquare.psth,NeuronStruct(cl).bigSquare.psth,NeuronStruct(cl).ffPulse.psth])])
    
    subplot(3,6,5)
    polar([NeuronStruct(cl).DriftGrat.slow.RadAngle,NeuronStruct(cl).DriftGrat.slow.RadAngle(1)],...
        [NeuronStruct(cl).DriftGrat.slow.SpikesNum,NeuronStruct(cl).DriftGrat.slow.SpikesNum(1)])
    hold on
    polar([NeuronStruct(cl).DriftGrat.slow.VectorAngle,NeuronStruct(cl).DriftGrat.slow.VectorAngle],...
        [0,NeuronStruct(cl).DriftGrat.slow.VectorSum],'r')
    title('slow drifting grating')
    
    subplot(3,6,6)
    polar([NeuronStruct(cl).DriftGrat.fast.RadAngle,NeuronStruct(cl).DriftGrat.fast.RadAngle(1)],...
        [NeuronStruct(cl).DriftGrat.fast.SpikesNum,NeuronStruct(cl).DriftGrat.fast.SpikesNum(1)])
    hold on
    polar([NeuronStruct(cl).DriftGrat.fast.VectorAngle,NeuronStruct(cl).DriftGrat.fast.VectorAngle],...
        [0,NeuronStruct(cl).DriftGrat.fast.VectorSum],'r')
    title('fast drifting grating')
    
    if saveFigFlag ;
         print(gcf, '-djpeg', [saveFigPath,'Cell',num2str(cl)])
    else
        disp(['cell ',num2str(cl),': ', NeuronStruct(cl).cellType]) ;
        InputTxt = input('cell type:','s') ;
        if ~isempty(InputTxt)
            NeuronStruct(cl).cellType = InputTxt ;
        end
    end
end
            
% assessing mosiacs of identified cell types
CellTypesList = {} ;
cnt=1 ;
for cl = 1:numCells ;
    if ~ismember(NeuronStruct(cl).cellType,CellTypesList) ; % if this cell type is not in the list
        CellTypesList{cnt} = NeuronStruct(cl).cellType ;
        cnt = cnt+1 ;
    end
end
    
for typ=1:length(CellTypesList) ; % for each cell type
    figure
    for cl = 1:numCells ; % for each cell
        if strcmp(CellTypesList{typ},NeuronStruct(cl).cellType) ; % if it is this cell type
            plot(NeuronStruct(cl).smallSquare.LikelyContour(1,:),NeuronStruct(cl).smallSquare.LikelyContour(2,:)) ; % plot is small contour
            hold on
        end
    end
    title(CellTypesList{typ})
    axis([0 400 0 400])
end
    
