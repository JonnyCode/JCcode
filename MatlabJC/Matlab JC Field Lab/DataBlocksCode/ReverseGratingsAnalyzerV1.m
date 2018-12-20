function ForIgor = ReverseGratingsAnalyzerV1(DataBlock, DB, Params)

% this function is adapted from
% 'KoReverseGratingsAnalyzerConcatPlusMappingMultiSetV1' to deal with no KO data

% JC 7/14/2017 

% parameters
Color_list = {'k','r','b','g','y','v'} ; % order of colors for each 
TrialTrigInterval = 10 ; % sec
InterTrialTime = 2 ; % sec
StimTransTime = .1 ; % (sec) skip cycles before this time after starting new stimulus
MicronPerPix = 4 ; % (um/pix)
EsemFlag = false; % calc psth esemble
FigureFlag = false ; % if you want to run figures
MseRelWashThreshold = 1 ; % (unitless) below this point response is "stable"
MseRelDrugDivWashThreshold = 1 ; % (unitless) above this point response has big drug effect 
LfSpThresh = 300 ; % (pix) spatial periods equal and above are termed "low frequency"

% load stimulus and get triggers for each part of concat data set

% load data 
dataRun = load_data(DataBlock(DB).RgPath{RgPathNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun, 'all') ; % load electrical images

% load stim
dataRun.names.stimulus_path = [DataBlock(DB).RgPath{RgPathNum}(1:end-15),'stimuli/s',DataBlock(DB).RgPath{RgPathNum}(end-1:end)] ;
dataRun = load_stim(dataRun,'user_defined_trigger_interval', TrialTrigInterval) ;

for t=1:length(dataRun.stimulus.trials) ; % for each trial
    trial_PHASE(t) = dataRun.stimulus.trials(t).SPATIAL_PHASE/dataRun.stimulus.trials(t).SPATIAL_PERIOD ;
end
params_PHASE = unique(trial_PHASE) ;

setTrialNum = length(dataRun.stimulus.trials) ; % number of trials for this stim set
setTrigNum = length(dataRun.triggers) ; % number of triggers in this data set

l_sp = length(dataRun.stimulus.params.SPATIAL_PERIOD) ;
l_tp = length(dataRun.stimulus.params.TEMPORAL_PERIOD) ;
l_phs = length(params_PHASE) ;

% stim trials with the same spatial and temporal periods
for sp=1:l_sp ; % for each spatial period
    for tp=1:l_tp ; % for each temporal period
        for phs=1:l_phs ; % for each contrast
            stim_groups_trials{sp,tp,phs} = [] ; % empty group
            for t=1:length(dataRun.stimulus.trials) ; % for each trial
                if dataRun.stimulus.trials(t).SPATIAL_PERIOD == dataRun.stimulus.params.SPATIAL_PERIOD(sp) ...
                        & dataRun.stimulus.trials(t).TEMPORAL_PERIOD == dataRun.stimulus.params.TEMPORAL_PERIOD(tp)...
                        & trial_PHASE(t) == params_PHASE(phs) ;
                    stim_groups_trials{sp,tp,phs} = [stim_groups_trials{sp,tp,phs}(:)',t] ; % this trial included
                end
            end
        end
    end
end

% trigger times for each (THIS ASSUMES the sp, tp,cntrst are all the
% same between data sets, and only the trigg times change)
for sp=1:l_sp ; % for each spatial period
    for tp=1:l_tp ; % for each temporal period
        for phs=1:l_phs ; % for each contrast
            stim_group_triggers_i{sp,tp,phs} = [] ;
            tempTrigs = dataRun.stimulus.triggers(stim_groups_trials{sp,tp,phs}) ; % the first trigger for each stimulus 
            for a = 1:length(tempTrigs) ; % for each first trigger
                stim_group_triggers_i{dset}{sp,tp,phs} = [stim_group_triggers_i{sp,tp,phs}(:); find(dataRun.triggers>=tempTrigs(a)+StimTransTime & ...
                    dataRun.triggers <= tempTrigs(a)+TrialTrigInterval-InterTrialTime)] ; % find the triggers within that set
            end
        end
    end
end
% stim params
SpatialPeriod = dataRun.stimulus.params.SPATIAL_PERIOD ;
TemporalPeriod = dataRun.stimulus.params.TEMPORAL_PERIOD ;
Phase = params_PHASE ;


Spatial_frequency = 1./(SpatialPeriod*MicronPerPix) ; % cycles/um

% mapping

l_cells = length(dataRun.spikes) ; % number of cells
% map using electrical images
cell_list_map = map_ei(dataRunMaster, dataRun) ;

% cells ids in slave for each UniqueCellType set in master data
for uc = 1:length(UniqueCellTypes) ;
    Tempi = get_cell_indices(dataRunMaster, UniqueCellTypes{uc}) ;
    cell_ids{uc} = cell2mat(cell_list_map(Tempi)) ;
    if ~isempty(cell_ids{uc}) ;
        cell_i{uc} = get_cell_indices(dataRun, cell_ids{uc}) ;
    else
        cell_i{uc} = [] ;
    end
end

% prep cells and get trigger times
for dset=1:length(DataBlock(DB).RgPath) ;
    for sp=1:l_sp ; % for each spatial period
        for phs=1:l_phs ; % for each contrast
            for tp=1:l_tp ; % for each temporal period
                stim_group_triggers{sp,tp,phs}= dataRun.triggers(stim_group_triggers_i{sp,tp,phs}) ; % trigger times
                psth{tp}{phs}{sp} = nans(l_cells,ceil(min(diff(stim_group_triggers{sp,tp,phs}'))*dataRun.sampling_rate)) ;
            end
        end
    end

    for tp=1:l_tp ; % for each temporal period
        for sp=1:l_sp ; % for each spatial period
            for cells = 1:l_cells ; % for each cell
                for phs=1:l_phs ; % for each phase
                    if ~isempty(stim_group_triggers{dset}{sp,tp,phs}) ;
                        [psthTemp,PsthTime] = get_smooth_psth(dataRun.spikes{cells},stim_group_triggers{dset}{sp,tp,phs}','stop',min(diff(stim_group_triggers{dset}{sp,tp,phs}'))) ;
                        psth{dset}{tp}{phs}{sp}(cells,:) = psthTemp ;

                        [powerspec_xvalues, mean_powerspec] = PowerSpectrumFinder(psth{dset}{tp}{phs}{sp}(cells,:),dataRun.sampling_rate) ;
                        psth_varF1{dset}{tp}{sp}(cells,phs) = mean_powerspec(2) ;
                        psth_varF2{dset}{tp}{sp}(cells,phs) = mean_powerspec(3) ;
                        
                        psth_var{dset}{tp}{sp}(cells,phs) = var(psth{dset}{tp}{phs}{sp}(cells,:)) ; 
                        psth_mean{dset}{tp}{sp}(cells,phs) = mean(psth{dset}{tp}{phs}{sp}(cells,:)) ;
                        psth_range{dset}{tp}{sp}(cells,phs) = range(psth{dset}{tp}{phs}{sp}(cells,:)) ;

                        psth_peak{dset}{tp}{sp}(cells,phs) = max(psth{dset}{tp}{phs}{sp}(cells,:)) ;
                        psth_peakDivMean{dset}{tp}{sp}(cells,phs) = psth_peak{dset}{tp}{sp}(cells,phs)/psth_mean{dset}{tp}{sp}(cells,phs) ;
                        psth_duty{dset}{tp}{sp}(cells,phs) = sum(psth{dset}{tp}{phs}{sp}(cells,:)>(psth_peak{dset}{tp}{sp}(cells,phs)/2))/length(psthTemp) ; % fraction of time above 50%  
                    end
                end
            end        
        end
    end
end

% averaged over phase
for dset=1:length(DataBlock(DB).RgPath) ;
    for tp=1:l_tp ; % for each temporal period
        for sp=1:l_sp ; % for each spatial period
            for cells = 1:l_cells ; % for each cell
                 psth_varF1_phsMean{dset}{tp}(cells,sp) = mean(psth_varF1{dset}{tp}{sp}(cells,:)) ;
                 psth_varF2_phsMean{dset}{tp}(cells,sp) = mean(psth_varF2{dset}{tp}{sp}(cells,:)) ;
                 psth_var_phsMean{dset}{tp}(cells,sp) = mean(psth_var{dset}{tp}{sp}(cells,:)) ;
            end
        end
    end
end

% cell list struct to mat
for a=1:length(cell_list_map) ;
    if isempty(cell_list_map{a}) ;
        cell_list_map_mat(a) = nan ;
    else
        cell_list_map_mat(a) = cell_list_map{a} ;
    end
end

% temporal receptive fields by cell type
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    if length(cell_ids{uc})>0 ; 
        for a=1:length(cell_ids{uc}) ; % for each cell of this type
            Mi = find(cell_list_map_mat==cell_ids{uc}(a)) ; % master cell index
            Mi = Mi(1) ; % some slaves seem to have two masters?
            if isempty(trf{Mi})
                disp(['empty trf: ',num2str(Mi)])
                trf_uc{uc}(a,:) = nan(size(trf_time)) ;
            else
                trf_uc{uc}(a,:) = trf{Mi}' ;
            end
        end
        trf_uc_mean(uc,:) = nanmean(trf_uc{uc},1) ;
    else
        trf_uc_mean(uc,:) = nan(size(trf_time)) ;
    end
end


% rename some datafields so can delete dataRun struct
dataRun_cell_ids = dataRun.cell_ids ; % so can delete dataRun
dataRunMaster_cell_ids = dataRunMaster.cell_ids ;
cell_type_list = get_cell_type_list(dataRunMaster) ;
coord_tform = coordinate_transform(dataRunMaster,'sta');
for c = 1:length(dataRunMaster.cell_ids) ; % for each with a BW mapped ref
    ctr{c} = dataRunMaster.stas.fits{c}.mean ;
    rad{c} = dataRunMaster.stas.fits{c}.sd ;
    angle{c} = dataRunMaster.stas.fits{c}.angle ;
    if ~isempty(cell_list_map{c}) ; % if there is a mapped cell
        slave_c(c) = get_cell_indices(dataRun,cell_list_map{c}) ;
    end
end

for uc = 1:length(UniqueCellTypes) ; % for each cell type
    Master_i{uc} = get_cell_indices(dataRunMaster,UniqueCellTypes{uc}) ;
end

clear dataRun dataRunMaster
save(['/Users/jcafaro/Desktop/Temp/Matfiles/KoRGratV1',num2str(DB)],'-v7.3')

% FIGURES
if FigureFlag ;
    
% psths of each cell in control, +psem, wash
cells=1 ;
while cells <=l_cells ; % for each cell
    figure(1)
    set(gcf,'name',[num2str(cells)])
    clf
    for dset=1:length(DataBlock(DB).RgPath) ;
        for sp=1:l_sp ; % for each spatial period 
            for phs=1:l_phs ; % for each contrast
                subplot(l_sp,l_phs,l_phs*(sp-1)+phs) ;
                plot(PsthTime,psth{dset}{1}{phs}{sp}(cells,:),'color',Color_list{dset})
                axis tight
                hold on
                title(['']) ;
            end
        end
    end
    
    figure(2)
    clf
    % master cell indicy
    masteri=[] ;
    for clm=1:length(cell_list_map) ;
        if cell_list_map{clm}==dataRun_cell_ids(cells) ;
            masteri=clm ;
        end
    end
    if ~isempty(masteri) && ~isempty(srf{masteri}) ; % if there is a mapped master cell that has a rf
        subplot(2,3,1) % spatial receptive field
        temp_rf = srf{masteri} ;
        norm_rf = norm_image(temp_rf);
        imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
        colormap(brewermap([],'RdBu'))
        caxis([0,1]) 
        set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
        title(num2str(dataRunMaster_cell_ids(masteri)))
        
         subplot(2,3,2) % temporal receptive field
        plot(fliplr(trf_time),trf{masteri}) ;
        title(cell_type_list{masteri})
    end
    for dset=1:length(DataBlock(DB).RgPath) ;
        subplot(2,3,3) % psth variance (spatial frequency)
        plot(log(Spatial_frequency),psth_varF1_phsMean{dset}{1}(cells,:),'color',Color_list{dset})
        hold on
        plot(log(Spatial_frequency),psth_varF2_phsMean{dset}{1}(cells,:),'color',Color_list{dset},'LineStyle','--')
        xlabel('log 2 spatial f')
        ylabel('psth variance')
        title(num2str(dataRun_cell_ids(cells)))

        subplot(2,3,4) % psth variance (spatial frequency)
        plot(log(Spatial_frequency),psth_var_phsMean{dset}{1}(cells,:),'color',Color_list{dset})
        xlabel('log 2 spatial f')
        ylabel('psth variance')
        title(num2str(dataRun_cell_ids(cells)))
        hold on

        subplot(2,3,5) % psth variance (spatial frequency)
        plot(log(Spatial_frequency),psth_varF1_phsMean{dset}{1}(cells,:),'color',Color_list{dset})
        xlabel('log 2 spatial f')
        ylabel('psth variance (F1)')
        title(num2str(dataRun_cell_ids(cells)))
        hold on

        subplot(2,3,6) % psth variance (spatial frequency)
        plot(log(Spatial_frequency),psth_varF2_phsMean{dset}{1}(cells,:),'color',Color_list{dset})
        xlabel('log 2 spatial f')
        ylabel('psth variance (F2)')
        title(num2str(dataRun_cell_ids(cells)))
        hold on
    end
    
    nxtval = input('next (cell num, 0=back, return=forward)') ;
    if isempty(nxtval) ;
        cells=cells+1 ;
    elseif nxtval == 0 ;
        cells=cells-1 ;
    elseif nxtval>0 ;
        cells=nxtval ;
    end
end

% psth [control, +PSEM, wash] organized by unique cell types for a single spatial freq and contrast
sp = 2 ; % spatial frequency default
phs = 1 ; % contrast default
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=[1,2,4] ; 
            plot(PsthTime,psth{dset}{1}{phs}{sp}(cell_i{uc}(cells),:),'color',Color_list{dset})
            hold on
        end
        title(num2str(cell_i{uc}(cells)))
    end
end


% contrast gain functions for each cell type
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=1:3 ; 
            plot(Spatial_frequency,psth_var_cgain{dset}{1}(cell_i{uc}(cells),:),'color',Color_list{dset})
            hold on
        end
        set(gca,'Xscale','log')
        title(num2str(cell_i{uc}(cells)))
    end
end
    
% spatial tuning functions averaged across contrasts for each cell type
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=1:length(DataBlock(DB).RgPath) ; 
            %plot(Spatial_frequency,psth_varF1_phsMean{dset}{1}(cell_i{uc}(cells),:),'color',Color_list{dset},'LineStyle','--')
            hold on 
            plot(Spatial_frequency,psth_varF2_phsMean{dset}{1}(cell_i{uc}(cells),:),'color',Color_list{dset})
        end
        set(gca,'Xscale','log')
        title(num2str(cell_i{uc}(cells)))
    end
end

% parameter changes as function of sp averaged across contrast
% organized by cell type
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        PlotColor = 'c' ;
        if MseRelWashLfMean(cell_i{uc}(cells))<MseRelWashThreshold ; % if this cell was stable
            PlotColor = 'b' ;
        end
        if MseRelDrug_DivMseRelWashLfMean(cell_i{uc}(cells))>MseRelDrugDivWashThreshold ; % if this cell had a drug change > wash change
            PlotColor = 'r' ;
        end
        if MseRelWashLfMean(cell_i{uc}(cells))<MseRelWashThreshold & MseRelDrug_DivMseRelWashLfMean(cell_i{uc}(cells))>MseRelDrugDivWashThreshold ; % if this cell was stable and had a bigger drug effect
            PlotColor = 'k' ;
        end
        
        subplot(5,1,1) 
        plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean(cell_i{uc}(cells),:), PlotColor)
        hold on
        
        subplot(5,1,2) 
        plot(Spatial_frequency,psth_var_cgain_DivCntrl{2}(cell_i{uc}(cells),:), PlotColor)
        hold on
        
        subplot(5,1,3) 
        plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean(cell_i{uc}(cells),:), PlotColor)
        hold on

        subplot(5,1,4) 
        plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean(cell_i{uc}(cells),:), PlotColor)
        hold on
    end
    
    subplot(5,1,1) 
    %plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean_AllCells(uc,:),'b','linewidth',2)
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'r:')
    set(gca,'Xscale','log')
    axis tight
    ylabel('var drug/cntrl')
    xlabel('Spatial Frequency (c/um)')

    subplot(5,1,2) 
    %plot(Spatial_frequency,psth_var_cgain_DrugDivCntrl_AllCells(uc,:),'b','linewidth',2)
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'r:')
    set(gca,'Xscale','log')
    axis tight
    ylabel('var gain drug/cntrl')
    xlabel('Spatial Frequency (c/um)')

    subplot(5,1,3) 
    %plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean_AllCells(uc,:),'b','linewidth',2)
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'r:')
    set(gca,'Xscale','log')
    axis tight
    ylabel('mean drug/cntrl')
    xlabel('Spatial Frequency (c/um)')

    subplot(5,1,4) 
    %plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean_AllCells(uc,:),'b','linewidth',2)
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'r:')
    set(gca,'Xscale','log')
    axis tight
    ylabel('duty drug/cntrl')
    xlabel('Spatial Frequency (c/um)')
end

% parameter changes as function of sp averaged across contrast
% comparing averages across cell types
figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    subplot(6,1,1) 
    plot(Spatial_frequency,nanmean(psth_var_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight

    subplot(6,1,2) 
    plot(Spatial_frequency,nanmean(psth_range_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight

    subplot(6,1,3) 
    plot(Spatial_frequency,nanmean(psth_mean_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight

    subplot(6,1,4) 
    plot(Spatial_frequency,nanmean(psth_duty_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight

    subplot(6,1,5) 
    plot(Spatial_frequency,nanmean(psth_peakDivMean_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight
end

% average stats for all cell types for 'stable cells' 
NumFigRows = 6 ;
figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    
    subplot(NumFigRows,length(UniqueCellTypes),uc) ; % trfs 
    for c = 1:length(Master_i{uc}) ; % for each cell
        [X,Y] = drawEllipse([ctr{Master_i{uc}(c)} rad{Master_i{uc}(c)} angle{Master_i{uc}(c)}]) ;
        if ~any(isnan([X,Y])) ;
            [X,Y] = tformfwd(coord_tform, X, Y) ;
            if isempty(cell_list_map{Master_i{uc}(c)})    
                plot(X,Y,'k')
            else
                plot(X,Y,'r')
            end
            hold on
        end
    end   
    title(['N=',num2str(length(cell_ids{uc})),'/',num2str(length(Master_i{uc}))])
    
    if cell_ids{uc}>0 ;
        Nsqrt = sqrt(length(cell_ids{uc})) ;

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*1+uc) ; % trfs
        plot(trf_time,trf_uc_mean(uc,:))
        axis tight
        title(UniqueCellTypes{uc}) 
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*2+uc) ; % variance (spatial frequecy)
        errorbar(Spatial_frequency,mean(psth_var_cMean{1}{1}(cell_i{uc}(cellsStable{uc}),:),1),std(psth_var_cMean{1}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{1}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,'color',Color_list{1})
        hold on
        errorbar(Spatial_frequency,mean(psth_var_cMean{2}{1}(cell_i{uc}(cellsStable{uc}),:),1),std(psth_var_cMean{2}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{2}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,'color',Color_list{2})
        errorbar(Spatial_frequency,mean(psth_var_cMean{3}{1}(cell_i{uc}(cellsStable{uc}),:),1),std(psth_var_cMean{3}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{3}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,'color',Color_list{3})
        set(gca,'Xscale','log')
        axis tight
        ylabel('var (sp/s)^2')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*3+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_var_WashdivCntrl_cMean_StableCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean_StableCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('var delta')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*4+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_mean_WashdivCntrl_cMean_StableCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean_StableCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('mean delta')
        xlabel('Sf (c/um)')
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*5+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_duty_WashdivCntrl_cMean_StableCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean_StableCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('duty delta')
        xlabel('Sf (c/um)')
    end
end
        
% cells with large Drug mse /wash mse 
NumFigRows = 6 ;
figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    
    subplot(NumFigRows,length(UniqueCellTypes),uc) ; % trfs
    for c = 1:length(Master_i{uc}) ; % for each cell
        [X,Y] = drawEllipse([ctr{Master_i{uc}(c)} rad{Master_i{uc}(c)} angle{Master_i{uc}(c)}]) ;
        if ~any(isnan([X,Y])) ;
            [X,Y] = tformfwd(coord_tform, X, Y) ;
            if isempty(cell_list_map{Master_i{uc}(c)})    
                plot(X,Y,'k')
            else
                plot(X,Y,'r')
            end
            hold on
        end
    end   
    title(['N=',num2str(length(cell_ids{uc})),'/',num2str(length(Master_i{uc}))])
    
    if cell_ids{uc}>0 ;
        Nsqrt = sqrt(length(cell_ids{uc})) ;

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*1+uc) ; % trfs
        plot(trf_time,trf_uc_mean(uc,:))
        axis tight
        title(UniqueCellTypes{uc}) 
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*2+uc) ; % variance (spatial frequecy)
        errorbar(Spatial_frequency,mean(psth_var_cMean{1}{1}(cell_i{uc}(cellsEffect{uc}),:),1),std(psth_var_cMean{1}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{1}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{1})
        hold on
        errorbar(Spatial_frequency,mean(psth_var_cMean{2}{1}(cell_i{uc}(cellsEffect{uc}),:),1),std(psth_var_cMean{2}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{2}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{2})
        errorbar(Spatial_frequency,mean(psth_var_cMean{3}{1}(cell_i{uc}(cellsEffect{uc}),:),1),std(psth_var_cMean{3}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{3}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{3})
        set(gca,'Xscale','log')
        axis tight
        ylabel('var (sp/s)^2')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*3+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_var_WashdivCntrl_cMean_EffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean_EffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('var delta')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*4+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_mean_WashdivCntrl_cMean_EffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean_EffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('mean delta')
        xlabel('Sf (c/um)')
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*5+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_duty_WashdivCntrl_cMean_EffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean_EffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('duty delta')
        xlabel('Sf (c/um)')
    end
end

% cells that are stable and have large drug/wash mse
NumFigRows = 6 ;
figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    
    subplot(NumFigRows,length(UniqueCellTypes),uc) ; % trfs
    for c = 1:length(Master_i{uc}) ; % for each cell
        [X,Y] = drawEllipse([ctr(Master_i{uc}(c)) rad(Master_i{uc}(c)) angle(Master_i{uc}(c))]) ;
        if ~any(isnan([X,Y])) ;
            [X,Y] = tformfwd(coord_tform, X, Y) ;
            if isempty(cell_list_map{Master_i(c)})    
                plot(X,Y,'k')
            else
                plot(X,Y,'r')
            end
            hold on
        end
    end   
    title(['N=',num2str(length(cell_ids{uc})),'/',num2str(length(Master_i{uc}))])
    
    if cell_ids{uc}>0 ;
        Nsqrt = sqrt(length(cell_ids{uc})) ;

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*1+uc) ; % trfs
        plot(trf_time,trf_uc_mean(uc,:))
        axis tight
        title(UniqueCellTypes{uc}) 
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*2+uc) ; % variance (spatial frequecy)
        errorbar(Spatial_frequency,mean(psth_var_cMean{1}{1}(cell_i{uc}(cellsStableEffect{uc}),:),1),std(psth_var_cMean{1}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{1}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{1})
        hold on
        errorbar(Spatial_frequency,mean(psth_var_cMean{2}{1}(cell_i{uc}(cellsStableEffect{uc}),:),1),std(psth_var_cMean{2}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{2}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{2})
        errorbar(Spatial_frequency,mean(psth_var_cMean{3}{1}(cell_i{uc}(cellsStableEffect{uc}),:),1),std(psth_var_cMean{3}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{3}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{3})
        set(gca,'Xscale','log')
        axis tight
        ylabel('var (sp/s)^2')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*3+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_var_WashdivCntrl_cMean_StableEffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean_StableEffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('var delta')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*4+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_mean_WashdivCntrl_cMean_StableEffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean_StableEffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('mean delta')
        xlabel('Sf (c/um)')
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*5+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_duty_WashdivCntrl_cMean_StableEffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean_StableEffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('duty delta')
        xlabel('Sf (c/um)')
    end
end
    
% For Igor export
ForIgor = struct ;

VecName = ['MseHistX','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,MseRelHistX) ; % mse bins

VecName = ['MseDrugOverWashCdf','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,MseRelDrug_DivMseRelWash_Hist_cumsum) ; % mse drug:wash ratio cumulative pdf

VecName = ['MseLfmDrugOverWashCdf','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,MseRelDrug_DivMseRelWashLfMean_Hist_cumsum) ; % mse drug:wash ratio cumulative pdf

VecName = ['MseLfmWashCdf','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,MseRelWashLfMean_Hist_cumsum) ; % mse drug:wash ratio cumulative pdf

VecName = ['MseLfmDrugCdf','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,MseRelDrugLfMean_Hist_cumsum) ; % mse drug:wash ratio cumulative pdf

exportStructToHDF5(ForIgor,['/Users/jcafaro/Desktop/Temp/Matfiles/KoSpatialV3Db',num2str(DB),'.h5'],'/');

else
ForIgor.Empty = [] ;
end



