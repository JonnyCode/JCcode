% DsMovSim3 Script to make plots

% datasets to load
DataSet{1} = 'RunId201812131545.mat' ;  % linear local
DataSet{2} = 'RunId201812131547.mat' ;  % linear random
DataSet{3} = 'RunId201812131549.mat' ;  % nonlinear local
DataSet{4} = 'RunId201812131551.mat' ;  % nonlinear random

SavePath = '/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/2018Paper/DsMovieSim3/mfile/' ;

figure % error as function of number of cells
ColorList = {'k:','r:','k','r'} ;
for Ds = 1:length(DataSet) ; % for each dataset
    
    % load
    load([SavePath,DataSet{Ds}],'NumOleCellsSet','ErrorMed_NumCellsLoop','ErrorMed_NumCellsLoop_Sem')

    % plot
    errorbar(NumOleCellsSet,ErrorMed_NumCellsLoop,ErrorMed_NumCellsLoop_Sem,ColorList{Ds})
    hold on

end
% save
SaveFigPath = '/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/2018Paper/DsMovieSim3Plots/figs/Plot1' ;
saveas(gcf,SaveFigPath)
print(gcf,'-dpdf',SaveFigPath)

