function ForIgor = DimFlashAnalysisWithCellTypes(DataBlock, DB, Params)

% This function assumes that cells are only typed in first data set.  Cells
% that are not mapped onto this are thrown out.

% JC 5/19/15
%DB=3; % TEMP

% parameters
mapEiFlag = DataBlock(DB).mapEiFlag ; % if datablocks should be mapped using map_ei
Color_list = {'k','r','g','y','b'} ; % order of colors for each 
interFlashIntVar = .005 ; % (sec) expected precision of trigger intervals
minNumFlashes = 5 ; % min number of flashes in a block 
NumStdaboveNoise = 2 ; % number of std above noise that is called signal

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

% load dim flash data
dataRun = load_data(DataBlock(DB).DfPath) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_params(dataRun) ;
dataRun = load_ei(dataRun, 'all') ;
                
% map dim flash onto bw data
if mapEiFlag ; % if using map ei cells
    % load electrical images
    dataRun = load_ei(dataRun, 'all') ;

    % map using electrical images
    cell_list_map = map_ei(dataRunMaster, dataRun, 'corr_threshold', 0.5) ;

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

% find sets of flash stimuli
PutativeSet_i = diff(dataRun.triggers)-DataBlock(DB).DfParams.interFlashInt < interFlashIntVar &...
    diff(dataRun.triggers)-DataBlock(DB).DfParams.interFlashInt > -interFlashIntVar ; % if time interval before next trigger is within range
st=1 ; 
trigger_set_i{1} = [] ;
for a=1:length(dataRun.triggers) ; % for each trigger
    trigger_set_i{st} = [trigger_set_i{st},a] ; % put it in a set   
    if a<length(dataRun.triggers) ; % if its not the last trigger
        if PutativeSet_i(a) ; % next trigger is has the right interval
            st = st ; % keep it in the set
        else
            st = st+1 ; % put it in a new set
            trigger_set_i{st} = [] ;
        end
    end
end

% get psth in each cell for each set of stimuli 
ds = 1 ; % identified trigger set
for ts=1:length(trigger_set_i) ; % for each possible trigger set  
    if length(trigger_set_i{ts})>minNumFlashes ; % if more than X flashes
        for uc = 1:length(UniqueCellTypes) ; % for each cell type
            if length(cell_ids{uc})>0 ; 
                for c = 1:length(cell_ids{uc}) ; % for each cell of that type
                    ci = get_cell_indices(dataRun, cell_ids{uc}(c)) ;
                    [psthTemp,psthTime] = get_smooth_psth(dataRun.spikes{ci}, dataRun.triggers(trigger_set_i{ts}), 'start', -2, 'stop', DataBlock(DB).DfParams.interFlashInt) ;
                    psth{uc}{ds}(c,:) = psthTemp ;
                    
                    preStimPnts = [1:find(psthTime==0)] ; % points before flash
                    preStimMean = mean(psthTemp(preStimPnts)) ; % mean before flash
                    preStimStd = std(psthTemp(preStimPnts)) ; % std before flash
                    
                    psth_delta{uc}{ds}(c,:) = psthTemp - preStimMean ; % response change
                    [v,i] = max(abs(psth_delta{uc}{ds}(c,:))) ; % largest response change
                    psth_delta_peak{uc}(ds,c) = psth_delta{uc}{ds}(c,i) ;
                    
                    RspPnts = find(abs(psth_delta{uc}{ds}(c,preStimPnts(end)+1:end))> preStimStd*NumStdaboveNoise) ; % "response points"
                    psth_absRspSum{uc}(ds,c) = sum(abs(psth_delta{uc}{ds}(c,preStimPnts(end)+RspPnts))) ; % sum of response points
                end 
                psth_mean{uc}(ds,:) = mean(psth{uc}{ds},1) ;
                psth_std{uc}(ds,:) = std(psth{uc}{ds},[],1) ;
            end
        end     
        ds = ds+1 ;
    end   
end

if length(psth{1}) > length(DataBlock(DB).DfParams.Ftime) ; % if the number of identified flash blocks is not as expected
    error('too many flash sets identified')
elseif length(psth{1}) < length(DataBlock(DB).DfParams.Ftime) ;
    error('too few flash sets identified')
end

% estimate the relative stim intensity assuming flashes are within linear response range
Irel = DataBlock(DB).DfParams.Ftime./10.^DataBlock(DB).DfParams.NDF ;
Irel_unique = unique(Irel) ;

% Organise by and average across cells with the same intensity levels 
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    if length(cell_ids{uc})>0 ; 
        psth_delta_peak_Irel{uc} = zeros(length(Irel_unique), length(cell_ids{uc})) ;
        psth_absRspSum_Irel{uc} = zeros(length(Irel_unique), length(cell_ids{uc})) ;

        for ds = 1:size(psth_mean{uc},1) ;
            [v,i] = intersect(Irel_unique,Irel(ds)) ; % index of unique intensity
            NumSame = sum(Irel_unique(i)==Irel) ; % number of same stim
            psth_delta_peak_Irel{uc}(i,:) = psth_delta_peak_Irel{uc}(i,:)+ psth_delta_peak{uc}(ds,:)/NumSame ;
            
            psth_absRspSum_Irel{uc}(i,:) = psth_absRspSum_Irel{uc}(i,:)+ psth_absRspSum{uc}(ds,:)/NumSame ;
        end
        psth_delta_peak_Irel_mean{uc} = mean(psth_delta_peak_Irel{uc},2) ;
        psth_delta_peak_Irel_sem{uc} = std(psth_delta_peak_Irel{uc},[],2)/sqrt(length(cell_ids{uc})) ;
        
        psth_absRspSum_Irel_mean{uc} = mean(psth_absRspSum_Irel{uc},2) ;
        psth_absRspSum_Irel_sem{uc} = std(psth_absRspSum_Irel{uc},[],2)/sqrt(length(cell_ids{uc})) ;
    end
end

% figures
Ftimes = unique(DataBlock(DB).DfParams.Ftime) ;
NDFs = unique(DataBlock(DB).DfParams.NDF) ;

for ds = 1:size(psth_mean{1},1) ;
    plot_i(ds) = (find(NDFs==DataBlock(DB).DfParams.NDF(ds))-1)*length(Ftimes) + find(Ftimes==DataBlock(DB).DfParams.Ftime(ds)) ; % subplot i
end

for uc = 1:length(UniqueCellTypes) ; % for each cell type
    if length(cell_ids{uc})>0 ; 
        figure
        for ds = 1:size(psth_mean{uc},1) ;
            subplot(length(NDFs),length(Ftimes),plot_i(ds))
            plot(psthTime,psth_mean{uc}(ds,:)) ;
            hold on
            %plot(psthTime,psth_mean{uc}(ds,:)+psth_std{uc}(ds,:),'b:') ;
            %plot(psthTime,psth_mean{uc}(ds,:)-psth_std{uc}(ds,:),'b:') ;
            text(.9,.9,['NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),',fdur ',num2str(DataBlock(DB).DfParams.Ftime(ds))],'units', 'norm') 
        end
        title(UniqueCellTypes{uc})
    end
end


figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    if length(cell_ids{uc})>0 ; 
        plot(log(Irel_unique),psth_delta_peak_Irel{uc})
        hold on
        plot(log(Irel_unique),psth_delta_peak_Irel_mean{uc},'k','linewidth',2)
    end
end

figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    if length(cell_ids{uc})>0 ; 
        plot(log(Irel_unique),psth_absRspSum_Irel{uc})
        hold on
        plot(log(Irel_unique),psth_absRspSum_Irel_mean{uc},'k','linewidth',2)
    end
end

% for Igor
ForIgor = struct() ;

VecName = ['psthTime','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psthTime) ; % psth time vector

for uc = 1:length(UniqueCellTypes) ; % for each cell type
    if length(cell_ids{uc})>0 ; 
        VecName = ['psth','Uc',num2str(uc),'Db',num2str(DB)] ;
        ForIgor = setfield(ForIgor, VecName, psth_mean{uc}') ; % psth vector
        
        VecName = ['psthPeak','Uc', num2str(uc),'Db',num2str(DB)] ;
        ForIgor = setfield(ForIgor, VecName, psth_delta_peak_Irel_mean{uc}) ; % peak vector
        
        VecName = ['psthPeakSEM','Uc', num2str(uc),'Db',num2str(DB)] ;
        ForIgor = setfield(ForIgor, VecName, psth_delta_peak_Irel_sem{uc}) ; % peak SEM vector
        
        VecName = ['N','Uc', num2str(uc),'Db',num2str(DB)] ;
        ForIgor = setfield(ForIgor, VecName, length(cell_ids{uc})) ; % peak vector
    end
end

VecName = ['Irel','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor, VecName, Irel_unique) ; % X axis I relative 

VecName = ['dsNDFandFtime','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor, VecName, [DataBlock(DB).DfParams.NDF;DataBlock(DB).DfParams.Ftime]) ; % psth vector






    