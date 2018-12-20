function ForIgor = SpontRateAnalysis(DataBlock, DB, Params)

% This function assumes that cells are only typed in first BW data set.  Cells
% that are not mapped onto this are thrown out.

% JC 6/19/15

% parameters
mapEiFlag = DataBlock(DB).mapEiFlag ; % if datablocks should be mapped using map_ei
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 
RateBinT_long = 1 ; % sec
RateBinT_short = .001 ; % sec
MaxLag = 1 ; % sec
Ac_Trange = [.05:.2] ; % (sec) look for max of AC between these values
numXrowsPlot = 5 ;

IsiDistX = [0:RateBinT_short:MaxLag] ;

% load Master data BW stimulus
dataRunMaster = load_data(DataBlock(DB).BwPath{1}) ;
dataRunMaster = load_neurons(dataRunMaster) ;
dataRunMaster = load_ei(dataRunMaster, 'all') ;
dataRunMaster = load_params(dataRunMaster) ;

for a = 1:length(dataRunMaster.cell_types) ;
    celltypes{a} = dataRunMaster.cell_types{a}.name ;
end

UniqueCellTypes = unique(celltypes) ;
if isempty(UniqueCellTypes{1}) ;
    UniqueCellTypes = UniqueCellTypes(2:end) ;
end

% load data
for SpN = 1:length(DataBlock(DB).SpontPath) ;
    dataRun = load_data(DataBlock(DB).SpontPath{SpN}) ;
    dataRun = load_neurons(dataRun) ;
    dataRun = load_params(dataRun) ;
    dataRun = load_ei(dataRun, 'all') ; 

    % map dim flash onto bw data
    if mapEiFlag ; % if using map ei cells
        % load electrical images
        dataRun = load_ei(dataRun, 'all') ;

        % map using electrical images
        cell_list_map = map_ei(dataRunMaster, dataRun, 'corr_threshold', 0.95) ;

        % cells ids in slave for each UniqueCellType set in master data
        for uc = 1:length(UniqueCellTypes) ;
            Tempi = get_cell_indices(dataRunMaster, UniqueCellTypes{uc}) ;
            cell_ids{uc} = cell2mat(cell_list_map(Tempi)) ;
        end
    else % if not using map ei
        for uc = 1:length(UniqueCellTypes) ;
            cell_ids{uc} = intersect(dataRun.cell_ids, get_cell_ids(dataRunMaster,UniqueCellTypes{uc})) ;
        end
    end
    
    % spike rate stats
    for uc = 1:length(UniqueCellTypes) ;
        if length(cell_ids{uc})>0 ;
            for c = 1:length(cell_ids{uc}) ; % for each cell of that type
                 ci = get_cell_indices(dataRun, cell_ids{uc}(c)) ;
                 SpikeRate{SpN}{uc}(c,:) = hist(dataRun.spikes{ci}, [0:RateBinT_long:dataRun.duration]) ;
                 SpikeRate_mean{SpN}{uc}(c) = mean(SpikeRate{SpN}{uc}(c,:)) ;
                 SpikeRate_std{SpN}{uc}(c) = mean(SpikeRate{SpN}{uc}(c,:)) ;   
            end
            SpikeRate_mean_mean{SpN}(uc) = mean(SpikeRate_mean{SpN}{uc}) ;
            SpikeRate_mean_sem{SpN}(uc) = std(SpikeRate_mean{SpN}{uc})/sqrt(length(cell_ids{uc})) ;
        end
    end
    
    % autocorrelation stats
    for uc = 1:length(UniqueCellTypes) ;
        if length(cell_ids{uc})>0 ;
            for c = 1:length(cell_ids{uc}) ; % for each cell of that type
                %IsiDist{SpN}{uc}(c,:) = hist(diff(dataRun.spikes{ci}),IsiDistX) ;         
                 st = zeros(1,dataRun.duration * dataRun.sampling_rate) ; % make a spike train
                 st(round(dataRun.spikes{ci} * dataRun.sampling_rate)) = 1 ; % populate spike train
                 Ac =  xcorr(st, MaxLag * dataRun.sampling_rate) ; % autocorrelation
                 Ac = smooth(smooth(Ac,RateBinT_short* dataRun.sampling_rate),RateBinT_short * dataRun.sampling_rate) ; % double smoothed Ac
                 Ac_smooth{SpN}{uc}(c,:) = Ac(MaxLag*dataRun.sampling_rate+1:end) ; % cut in half
                 Ac_ratio{SpN}{uc}(c) = max(Ac_smooth{SpN}{uc}(c,Ac_Trange*dataRun.sampling_rate))/Ac_smooth{SpN}{uc}(c,1) ; % peak ratio
            end
        end
    end 
 end

% figures

figure
for SpN = 1:length(DataBlock(DB).SpontPath) ;
    for uc = 1:length(UniqueCellTypes) ;
        if length(cell_ids{uc})>0 ;
            for c = 1:length(cell_ids{uc}) ; % for each cell of that type
                plot(SpikeRate{SpN}{uc}(c,:))
                pause
            end
        end
    end
end



figure 
for SpN = 1:length(DataBlock(DB).SpontPath) ;
    for uc = 1:length(UniqueCellTypes) ;
        subplot(numXrowsPlot,ceil(length(UniqueCellTypes)/numXrowsPlot),uc) 
        errorbar([1:length(cell_ids{uc})],SpikeRate_mean{SpN}{uc}, SpikeRate_std{SpN}{uc},SpikeRate_std{SpN}{uc},['*',Color_list{SpN}]) ;
        hold on
        errorbar(0, SpikeRate_mean_mean{SpN}(uc), SpikeRate_mean_sem{SpN}(uc), SpikeRate_mean_sem{SpN}(uc),['o',Color_list{SpN}])
        title(UniqueCellTypes{uc})
    end
end

% forIgor
ForIgor = struct() ;

for SpN = 1:length(DataBlock(DB).SpontPath) ;
    for uc = 1:length(UniqueCellTypes) ;  % for each cell type
        
        VecName = ['SpontRateMeans','Uc',num2str(uc),'Spn',num2str(SpN),'Db',num2str(DB)] ;
        ForIgor = setfield(ForIgor,VecName,SpikeRate_mean{SpN}{uc}) ; % time course

        VecName = ['SpontRateStd','Uc',num2str(uc),'Spn',num2str(SpN),'Db',num2str(DB)] ;
        ForIgor = setfield(ForIgor,VecName,SpikeRate_std{SpN}{uc}) ; % time course
        
        VecName = ['SpontRateMeanAllCells','Uc',num2str(uc),'Spn',num2str(SpN),'Db',num2str(DB)] ;
        ForIgor = setfield(ForIgor,VecName,SpikeRate_mean_mean{SpN}(uc)) ; % time course

        VecName = ['SpontRateSemAllCells','Uc',num2str(uc),'Spn',num2str(SpN),'Db',num2str(DB)] ;
        ForIgor = setfield(ForIgor,VecName,SpikeRate_mean_sem{SpN}(uc)) ; % time course
    end
end





