function ForIgor = MovieSpikeParamFinder(DataBlock, DB, Params)

% this function will calculate a set of spike parameters across a set of
% movies mapped from a BWN. Adapted from 'MovieRasterPlots.m'.

% JC 10/17/2016 

% parameters 
psthBin = .1 ; %(s)

% binary white data run
% load Master data
dataRunMaster = load_data(DataBlock(DB).BwPath) ;
dataRunMaster = load_neurons(dataRunMaster) ;
dataRunMaster = load_ei(dataRunMaster, 'all') ;
dataRunMaster = load_params(dataRunMaster,'cell_type_depth', 7) ;

dbl = length(DataBlock(DB).MovieRepPath) ;
numCells = length(dataRunMaster.spikes) ;

% load and map each movie block data set
for db=1:dbl ; % for each movie data block
    
    dataRun{db} = load_data(DataBlock(DB).MovieRepPath{db}) ;
    dataRun{db} = load_neurons(dataRun{db}) ;
    dataRun{db} = load_ei(dataRun{db}, 'all') ;
    
    cell_list_map{db} = map_ei(dataRunMaster, dataRun{db}) ;
    
    trigNum(db) = ceil(DataBlock(DB).MovieRepFrameNum(db)/100) ; % number of triggers per repeat
    trigs{db} = dataRun{db}.triggers(1:trigNum:end) ; % 
    repeatNum(db) = floor(length(dataRun{db}.triggers)/trigNum(db)); % show how many repeats this is
end

% params
stim_time =min(DataBlock(DB).MovieRepFrameNum)/60 ;

for cells = 1:numCells ; % for each cell
    for db=1:dbl ; % for each movie data block
        cell_id = cell_list_map{db}(cells) ;
        if ~isempty(cell_id{1}) ;
            cell_i = get_cell_indices(dataRun{db},cell_id{1}) ;
            for tr=1:length(trigs{db}) ;
                spikeT = dataRun{db}.spikes{cell_i}-trigs{db}(tr) ;
                psth = hist(spikeT(spikeT>=0 & spikeT<=stim_time),[0:psthBin:stim_time]) ;
                psth = psth/psthBin ;
                psth_mean{db}(cells,tr) = mean(psth) ;
                psth_max{db}(cells,tr) = max(psth) ;
                psth_sparsityIndex{db}(cells,tr) = mean(psth)^2/mean(psth.^2) ; %(0-1 lower is more sparse; Ahmed et al 2009, "hippocampal rate code..")
            end
            psth_mean_TrialMean(db,cells) = mean(psth_mean{db}(cells,:)) ;
            psth_mean_TrialStd(db,cells) = std(psth_mean{db}(cells,:)) ;
            psth_max_TrialMean(db,cells) = mean(psth_max{db}(cells,:)) ;
            psth_max_TrialStd(db,cells) = mean(psth_max{db}(cells,:)) ;
            psth_sparsityIndex_TrialMean(db,cells) = mean(psth_sparsityIndex{db}(cells,:)) ;
            psth_sparsityIndex_TrialStd(db,cells) = mean(psth_sparsityIndex{db}(cells,:)) ;
        else
            for tr=1:length(trigs{db}) ;
                psth_mean{db}(cells,tr) = nan ;
                psth_max{db}(cells,tr) = nan ;
                psth_sparsityIndex{db}(cells,tr) = nan ;
            end
            psth_mean_TrialMean(db,cells) = nan;
            psth_mean_TrialStd(db,cells) = nan ;
            psth_max_TrialMean(db,cells) = nan ;
            psth_max_TrialStd(db,cells) = nan ;
            psth_sparsityIndex_TrialMean(db,cells) = nan ;
            psth_sparsityIndex_TrialStd(db,cells) = nan ;
        end      
    end
end

% organize by cell types 
for a = 1:length(dataRunMaster.cell_types) ;
    celltypes{a} = dataRunMaster.cell_types{a}.name ;
end

UniqueCellTypes = unique(celltypes) ;
if isempty(UniqueCellTypes{1}) ;
    UniqueCellTypes = UniqueCellTypes(2:end) ;
end

% cells ids in slave for each UniqueCellType set in master data
for uc = 1:length(UniqueCellTypes) ;
    Cell_i{uc} = get_cell_indices(dataRunMaster, UniqueCellTypes{uc}) ; % cell indicy of master

    psth_mean_TrialMean_CellTypeOrg{uc} = psth_mean_TrialMean(:,Cell_i{uc}) ;
    psth_max_TrialMean_CellTypeOrg{uc} = psth_max_TrialMean(:,Cell_i{uc}) ;
    psth_sparsityIndex_TrialMean_CellTypeOrg{uc} = psth_sparsityIndex_TrialMean(:,Cell_i{uc}) ;
end

% figures

figure % psth mean for each cell
set(gcf,'name','mean rate')
for db=1:dbl ; % for each movie data block
    subplot(dbl+1,3,3*db-2)
    errorbar([1:numCells],psth_mean_TrialMean(db,:),psth_mean_TrialStd(db,:),'*')
    axis tight
    ylabel(['mov',num2str(db)])
end
   
subplot(dbl+1,3,3*(dbl+1)-2)
errorbar([1:numCells],nanmean(psth_mean_TrialMean,1),nanstd(psth_mean_TrialMean,[],1),'*')
ylabel('average')
xlabel('cell i')
axis tight

for db=1:dbl ; % for each movie data block
    for uc = 1:length(UniqueCellTypes) ; % for each cell type
        subplot(dbl+1,3,3*db-1)
        plot(ones(1,length(Cell_i{uc}))*uc,psth_mean_TrialMean(db,Cell_i{uc}),'*')
        hold on
        axis tight
    end
end

for uc = 1:length(UniqueCellTypes) ; % for each cell type
    subplot(dbl+1,3,3*(dbl+1)-1)
    errorbar(uc,nanmean(psth_mean_TrialMean_CellTypeOrg{uc}(:)),nanstd(psth_mean_TrialMean_CellTypeOrg{uc}(:)),'*')
    hold on
    axis tight
end
xlabel('cell type')

for uc = 1:length(UniqueCellTypes) ; % for each cell type
    subplot(dbl+1,3,[3:3:3*dbl+1])
    plot(nanmean(psth_mean_TrialMean_CellTypeOrg{uc},2),[1:dbl],'*-')
    hold on
    axis tight
end
set(gca,'yDir','reverse')
ylabel('movie')
xlabel('mean rate')
    


figure % psth max for each cell
set(gcf,'name','max rate')
for db=1:dbl ; % for each movie data block
    subplot(dbl+1,3,3*db-2)
    errorbar([1:numCells],psth_max_TrialMean(db,:),psth_max_TrialStd(db,:),'*')
    axis tight
    ylabel(['mov',num2str(db)])
end
   
subplot(dbl+1,3,3*(dbl+1)-2)
errorbar([1:numCells],nanmean(psth_max_TrialMean,1),nanstd(psth_max_TrialMean,[],1),'*')
ylabel('average')
xlabel('cell i')
axis tight

for db=1:dbl ; % for each movie data block
    for uc = 1:length(UniqueCellTypes) ; % for each cell type
        subplot(dbl+1,3,3*db-1)
        plot(ones(1,length(Cell_i{uc}))*uc,psth_max_TrialMean(db,Cell_i{uc}),'*')
        hold on
        axis tight
    end
end

for uc = 1:length(UniqueCellTypes) ; % for each cell type
    subplot(dbl+1,3,3*(dbl+1)-1)
    errorbar(uc,nanmean(psth_max_TrialMean_CellTypeOrg{uc}(:)),nanstd(psth_max_TrialMean_CellTypeOrg{uc}(:)),'*')
    hold on
    axis tight
end
xlabel('cell type')

for uc = 1:length(UniqueCellTypes) ; % for each cell type
    subplot(dbl+1,3,[3:3:3*dbl+1])
    plot(nanmean(psth_max_TrialMean_CellTypeOrg{uc},2),[1:dbl],'*-')
    hold on
    axis tight
end
set(gca,'yDir','reverse')
ylabel('movie')
xlabel('max rate')
    

figure % psth sparsityIndex for each cell
set(gcf,'name','sparsity index')
for db=1:dbl ; % for each movie data block
    subplot(dbl+1,3,3*db-2)
    errorbar([1:numCells],psth_sparsityIndex_TrialMean(db,:),psth_sparsityIndex_TrialStd(db,:),'*')
    axis tight
    ylabel(['mov',num2str(db)])
end
   
subplot(dbl+1,3,3*(dbl+1)-2)
errorbar([1:numCells],nanmean(psth_sparsityIndex_TrialMean,1),nanstd(psth_sparsityIndex_TrialMean,[],1),'*')
ylabel('average')
xlabel('cell i')
axis tight

for db=1:dbl ; % for each movie data block
    for uc = 1:length(UniqueCellTypes) ; % for each cell type
        subplot(dbl+1,3,3*db-1)
        plot(ones(1,length(Cell_i{uc}))*uc,psth_sparsityIndex_TrialMean(db,Cell_i{uc}),'*')
        hold on
        axis tight
    end
end

for uc = 1:length(UniqueCellTypes) ; % for each cell type
    subplot(dbl+1,3,3*(dbl+1)-1)
    errorbar(uc,nanmean(psth_sparsityIndex_TrialMean_CellTypeOrg{uc}(:)),nanstd(psth_sparsityIndex_TrialMean_CellTypeOrg{uc}(:)),'*')
    hold on
    axis tight
end
xlabel('cell type')

for uc = 1:length(UniqueCellTypes) ; % for each cell type
    subplot(dbl+1,3,[3:3:3*dbl+1])
    plot(nanmean(psth_sparsityIndex_TrialMean_CellTypeOrg{uc},2),[1:dbl],'*-')
    hold on
    axis tight
end
set(gca,'yDir','reverse')
ylabel('movie')
xlabel('max rate')

