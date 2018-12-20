function ForIgor = DfAnalyzerV3WithBWmapping(DataBlock, DB, Params)

% adapted from DfAnalyzerV3 to include selection for particular cell types in BW
% also trimmed fat and made psth from all trial blocks

% JC 1/18/17

% parameters
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 
interFlashIntVar = .005 ; % (sec) expected precision of trigger intervals
minNumFlashes = 5 ; % min number of flashes in a block 
psth_bin_size = 0.1 ; % (sec) 
psth_start_time = -2 ; %(sec)
StdAbove = 2 ; % number of std above mean
LED_settings_calibrations_path = '/Volumes/lab/Experiments/Calibration/LED_settings_calibrations' ; % location of LED calibration mat file
LoadCalibrationFlag = true ;
PdfHistFlag = true ; % make histograms in pdfs if true
mapEiFlag = true ;
BwBlock = DataBlock(DB).DfBwBlock ; % block of BW to use as map
ExampleCellType = 'OnT' ; % cell type name that will be searched more specifically and exported

% estimate the relative stim intensity assuming flashes are within linear response range
if isfield(DataBlock(DB).DfParams,'Setting') ;
    load(LED_settings_calibrations_path) 
    for a=1:length(DataBlock(DB).DfParams.Ftime) ;
        ci = find(LED_settings_calibrations.FlashDuration==DataBlock(DB).DfParams.Ftime(a)) ;
        ri = find(round(LED_settings_calibrations.setting*10)==round(DataBlock(DB).DfParams.Setting(a)*10)) ;
        Irel(a) = LED_settings_calibrations.NormPower(ri,ci)*DataBlock(DB).DfParams.Ftime(a)./10.^DataBlock(DB).DfParams.NDF(a) ;
    end
elseif LoadCalibrationFlag ;
    load /Volumes/lab/Experiments/Calibration/NdfCalibration
    Irel = (DataBlock(DB).DfParams.Ftime/1000).*NdfCalibration(2,DataBlock(DB).DfParams.NDF+1) ;
else
    Irel = DataBlock(DB).DfParams.Ftime./10.^DataBlock(DB).DfParams.NDF ;
end
Irel_unique = unique(Irel) ;

% load data for dim flashes
dataRun = load_data(DataBlock(DB).DfPath) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun, 'all') ;
%dataRun = load_params(dataRun) ;

% params
NumCells = length(dataRun.spikes) ;

% find sets of flash stimuli
st=1 ; 
trigger_set_i{1} = [] ;
for a=1:length(dataRun.triggers) ; % for each trigger
    trigger_set_i{st} = [trigger_set_i{st},a] ; % put it in a set   
    if a<length(dataRun.triggers) ; % if its not the last trigger
        if sum(abs(dataRun.triggers(a+1)-dataRun.triggers(a)-DataBlock(DB).DfParams.interFlashInt)<interFlashIntVar)>0 ; % next trigger is has the right interval
            st = st ; % keep it in the set
        else
            st = st+1 ; % put it in a new set
            trigger_set_i{st} = [] ;
        end
    end
end

% select trigger sets that have enough trials to be counted
ds= 0 ;
for ts=1:length(trigger_set_i) ; % for each trigger set 
    if length(trigger_set_i{ts})>minNumFlashes ; % if more than X flashes
        ds=ds+1 ;
        select_trigger_set_i{ds} = trigger_set_i{ts} ;
    end
end

if length(select_trigger_set_i) ~= length(DataBlock(DB).DfParams.Ftime) ; % if the number of identified flash blocks is not as expected
    disp('Params or data length mismatch')
end

% group trigger sets with the same stimuli and order by stim intensity
for ds = 1:length(Irel_unique) ; % for each stimulus
    full_trigger_set_i{ds} = [] ;
    Ireli = find(Irel == Irel_unique(ds)) ; % find the 
    for ts=1:length(Ireli) ;
        full_trigger_set_i{ds} = [full_trigger_set_i{ds},select_trigger_set_i{Ireli(ts)}] ;
    end
end

UniqueStimNum = length(Irel_unique) ; % number of different flashes

% get psth in each cell for each set of stimuli 
for ts=1:UniqueStimNum ; % for each trigger set  
    for c = 1:NumCells ; % for each cell
        [psthTemp,psthTime] = get_smooth_psth(dataRun.spikes{c},...
            dataRun.triggers(full_trigger_set_i{ts}), 'start', psth_start_time, ...
            'stop', DataBlock(DB).DfParams.interFlashInt, 'bin_size', psth_bin_size) ; 
        psth{ts}(c,:) = psthTemp ;

        StimPnt = find(psthTime==0) ; % flash point
        preStimPnts = [1:StimPnt] ; % points before flash
        postStimPnts = [StimPnt+1:length(psthTime)] ; % points after flash

        preStimMean = mean(psthTemp(preStimPnts)) ; % mean before flash
        preStimStd = std(psthTemp(preStimPnts)) ; % std before flash
        preStimRange = range(psthTemp(preStimPnts)) ; % max before flash

        psth_delta{ts}(c,:) = psthTemp - preStimMean ; % response change
        [v,i] = max(abs(psth_delta{ts}(c,postStimPnts))) ; % largest response change after flash
        psth_delta_peak(ts,c) = psth_delta{ts}(c,postStimPnts(i)) ;

        psth_delta_peak_abs(ts,c) = abs(psth_delta_peak(ts,c)) ;
        psth_delta_peak_abs_Snr(ts,c) = psth_delta_peak_abs(ts,c)/preStimStd ;

        psth_delta_peak_time(ts,c) = psthTime(postStimPnts(i)) ; % time of peak after stim onset

        i1 = find(abs(psth_delta{ts}(c,1:postStimPnts(i)))<(v/2),1,'last') ; % first point
        i2 = find(abs(psth_delta{ts}(c,postStimPnts(i):end))<v/2,1,'first')+postStimPnts(i) ; % last point 
        if ~isempty(i1 + i2) ;
            psth_delta_halfWidth(ts,c) = psthTime(i2) - psthTime(i1) ; % half width
        end

        for trls = 1:length(full_trigger_set_i{ts}) ;
            trigTime = dataRun.triggers(full_trigger_set_i{ts}(trls))+psthTime(postStimPnts(i)) ;

            psth_peak_trls{ts}(c,trls) = sum(((trigTime-psth_bin_size/2)<dataRun.spikes{c}).*(dataRun.spikes{c}<(trigTime+psth_bin_size/2)))/psth_bin_size ;    
            Bg_rate_trls{ts}(c,trls) = sum(((dataRun.triggers(full_trigger_set_i{ts}(trls))-1)<dataRun.spikes{c}).*(dataRun.spikes{c}<dataRun.triggers(full_trigger_set_i{ts}(trls)))) ;
            psth_peak_delta_trls{ts}(c,trls) = psth_peak_trls{ts}(c,trls)- Bg_rate_trls{ts}(c,trls) ;
        end

        psth_peak_trls_mean(ts,c) = mean(psth_peak_trls{ts}(c,:),2) ;
        Bg_rate_trls_mean(ts,c) = mean(Bg_rate_trls{ts}(c,:),2) ;

        psth_peak_delta_trls_mean(ts,c) = mean(psth_peak_delta_trls{ts}(c,:),2) ;
        psth_peak_delta_trls_std(ts,c) = std(psth_peak_delta_trls{ts}(c,:),[],2) ;
    end
    psth_mean(ts,:) = mean(psth{ts}) ; % average across all cells
    psth_std(ts,:) = std(psth{ts}) ;

    psth_delta_mean(ts,:) = mean(psth_delta{ts}) ; 
end

Bg_rate_trls_mean_meanAcrossStim = mean(Bg_rate_trls_mean,1) ; % mean across stim (ie. 1 bg rate for each cell)

% Thresholds
NoThreshNum = 0 ; 
for c = 1:NumCells ; % for each cell
    ti = find(abs(psth_peak_delta_trls_mean(:,c))>(psth_peak_delta_trls_std(:,c)*StdAbove),1,'first') ;
    if ~isempty(ti) ;
        Threshold(c) = log10(Irel_unique(ti)) ;
    else
        Threshold(c) = nan ;
        NoThreshNum = NoThreshNum+1 ;
    end
end
FracNoThresh = NoThreshNum/NumCells ; % fraction of cells with no detected threshold

% peak rates
last_peak = psth_peak_delta_trls_mean(end,:) ; % max of tuning curve 

% half width and peak time of the response to the max stim
last_peak_time = psth_delta_peak_time(end,:) ;

last_halfWidth = psth_delta_halfWidth(end,:) ;

% histograms
peak_HistBins = [-200:5:200] ;
for ds = 1:size(psth_mean,1) ;
    peak_HistValues(ds,:) = cumsum(hist(psth_delta_peak(ds,:),peak_HistBins)) ; 
end

Threshold_HistBins = [-3:.1:3] ;
Threshold_HistValues= cumsum(hist(Threshold,Threshold_HistBins)) ; 

last_peak_HistBins = peak_HistBins ;
last_peak_HistValues = cumsum(hist(last_peak,last_peak_HistBins)) ; 

last_peak_time_HistBins = 0:.01:1 ;  
last_peak_time_HistValues = cumsum(hist(last_peak_time,last_peak_time_HistBins)) ;

last_halfWidth_HistBins = 0:.01:1 ;  
last_halfWidth_HistValues = cumsum(hist(last_halfWidth,last_halfWidth_HistBins)) ;

Bg_rate_HistBins = [0:1:100] ;
Bg_rate_HistValues = cumsum(hist(Bg_rate_trls_mean_meanAcrossStim,Bg_rate_HistBins)) ;

if PdfHistFlag ;
    Threshold_HistValues = Threshold_HistValues/Threshold_HistValues(end) ;
    Bg_rate_HistValues = Bg_rate_HistValues/Bg_rate_HistValues(end) ;
end
    
% map from BW
% load Master data BW stimulus
marks_params.thresh = 3 ;
dataRunMaster = load_data(DataBlock(DB).BwPath{BwBlock}) ;
dataRunMaster = load_neurons(dataRunMaster) ;
dataRunMaster = load_ei(dataRunMaster, 'all') ;
dataRunMaster = load_params(dataRunMaster) ;
% dataRunMaster = load_sta(dataRunMaster) ;
% dataRunMaster = get_sta_summaries(dataRunMaster, 'all','marks_params', marks_params) ;
% dataRunMaster = get_sta_fits_from_vision(dataRunMaster) ;
% dataRunMaster = get_autocorrelations(dataRunMaster, 'all') ;

for a = 1:length(dataRunMaster.cell_types) ;
    celltypes{a} = dataRunMaster.cell_types{a}.name ;
end

UniqueCellTypes = unique(celltypes) ;
if isempty(UniqueCellTypes{1}) ;
    UniqueCellTypes = UniqueCellTypes(2:end) ;
end

% map dim flash onto bw data
if mapEiFlag ; % if using map ei cells

    % map using electrical images
    cell_list_map = map_ei(dataRunMaster, dataRun) ;

    % cells ids in slave for each UniqueCellType set in master data
    for uc = 1:length(UniqueCellTypes) ;
        Masteri{uc} = get_cell_indices(dataRunMaster, UniqueCellTypes{uc}) ;
        cell_ids{uc} = cell2mat(cell_list_map(Masteri{uc})) ;
    end
else % if not using map ei
    for uc = 1:length(UniqueCellTypes) ;
        cell_ids{uc} = intersect(dataRun.cell_ids, get_cell_ids(dataRunMaster,UniqueCellTypes{uc})) ;
    end
end

% find histograms for each cell type
for uc = 1:length(UniqueCellTypes) ; 
    cell_i{uc} = get_cell_indices(dataRun,cell_ids{uc}) ;
    Threshold_byCellType_HistValues{uc} = cumsum(hist(Threshold(cell_i{uc}),Threshold_HistBins)) ;
    FractNoThresh_byCellType{uc} = sum(isnan(Threshold(cell_i{uc})))/length(cell_i{uc}) ;
    Bg_rate_byCellType_HistValues{uc} = cumsum(hist(Bg_rate_trls_mean_meanAcrossStim(cell_i{uc}),Bg_rate_HistBins)) ;
    last_peak_byCellType_HistValues{uc} = cumsum(hist(last_peak(cell_i{uc}),last_peak_HistBins)) ; 
    last_peak_time_byCellType_HistValues{uc} = cumsum(hist(last_peak_time(cell_i{uc}),last_peak_time_HistBins)) ;
    last_halfWidth_byCellType_HistValues{uc} = cumsum(hist(last_halfWidth(cell_i{uc}),last_halfWidth_HistBins)) ;
end
 
% example cell type index
ExampleCellTypei = nan ;
for uc = 1:length(UniqueCellTypes) ;
    if strcmp(ExampleCellType,UniqueCellTypes{uc})
        ExampleCellTypei = uc ;
    end
end
    
% figures
% 
% 
% % classification figures
% coord_tform = coordinate_transform(dataRunMaster,'sta');
% for c = 1:length(dataRunMaster.cell_ids) ; % for each with a BW mapped ref
%     ctr{c} = dataRunMaster.stas.fits{c}.mean ;
%     rad{c} = dataRunMaster.stas.fits{c}.sd ;
%     angle{c} = dataRunMaster.stas.fits{c}.angle ;
% end
% 
% dataRunMaster.stimulus.monitor_refresh = 60.35 ; % the monitor_refresh rate is not accurate
% time_course_time = -[0:dataRunMaster.stas.depth-1]*dataRunMaster.stimulus.interval/dataRunMaster.stimulus.monitor_refresh ;
% 
% figure % mosaic of example cell type
% for c = 1:length(Masteri{ExampleCellTypei}) ; % for each cell of the example type
%     [X,Y] = drawEllipse([ctr{Masteri{ExampleCellTypei}(c)}, ...
%         rad{Masteri{ExampleCellTypei}(c)}, angle{Masteri{ExampleCellTypei}(c)}]) ;
%     if ~any(isnan([X,Y])) ;
%         [X,Y] = tformfwd(coord_tform, X, Y) ;
%         if isempty(cell_list_map{Masteri{ExampleCellTypei}(c)})    
%             plot(X,Y,'k')
%         else
%             plot(X,Y,'r')
%         end
%         hold on
%     end
% end   
%     
% figure % time courses of example cell and all other cells
% for c = 1:length(dataRunMaster.cell_ids) ; % for each with a BW mapped ref
%     plot(fliplr(time_course_time), dataRunMaster.stas.time_courses{c}/...
%         norm(dataRunMaster.stas.time_courses{c}),'color',[.7,.7,.7])
%     hold on
% end
% Mi = get_cell_indices(dataRunMaster,ExampleCellType) ;
% for c = 1:length(Mi) ; % for each cell of the example cell type
%      plot(fliplr(time_course_time), dataRunMaster.stas.time_courses{Mi(c)}/...
%          norm(dataRunMaster.stas.time_courses{Mi(c)}),'k')
% end
%     
% figure % autocorrelations of example cell and all other cells
% for c = 1:length(dataRunMaster.cell_ids) ; % for each with a BW mapped ref
%     plot(dataRunMaster.autocorrelation{c}.bins, ...
%         dataRunMaster.autocorrelation{c}.probabilities/...
%         norm(dataRunMaster.autocorrelation{c}.probabilities),'color',[.7,.7,.7])
%     hold on
% end
% Mi = get_cell_indices(dataRunMaster,ExampleCellType) ;
% for c = 1:length(Mi) ; % for each cell of the example cell type
%      plot(dataRunMaster.autocorrelation{Mi(c)}.bins, ...
%          dataRunMaster.autocorrelation{Mi(c)}.probabilities/...
%          norm(dataRunMaster.autocorrelation{Mi(c)}.probabilities),'k')
% end
% 
figure % all cell histograms
plot(Threshold_HistBins,Threshold_HistValues) ;
hold on
plot(Threshold_HistBins(1),FracNoThresh,'r*')
ylabel('fraction of cells')
xlabel('Rh*/rod')

% figure % by cell type histograms
% numColumns = 3 ;
% for uc = 1:length(UniqueCellTypes) ;
%     subplot(ceil(length(UniqueCellTypes)/numColumns),numColumns,uc)
%     plot(Threshold_HistBins(1:end-1),Threshold_byCellType_HistValues{uc}(1:end-1)) ;
%     hold on
%     plot(Threshold_HistBins(1),Threshold_byCellType_HistValues{uc}(end),'r*')
%     title(UniqueCellTypes{uc})
% end
% 
% figure % by cell type histograms
% numColumns = 3 ;
% for uc = 1:length(UniqueCellTypes) ;
%     subplot(ceil(length(UniqueCellTypes)/numColumns),numColumns,uc)
%     plot(last_halfWidth_HistBins,last_halfWidth_byCellType_HistValues{uc}) ;
%     title(UniqueCellTypes{uc})
% end
% 
% figure % by cell type histograms
% numColumns = 3 ;
% for uc = 1:length(UniqueCellTypes) ;
%     subplot(ceil(length(UniqueCellTypes)/numColumns),numColumns,uc)
%     plot(last_peak_HistBins,last_peak_byCellType_HistValues{uc}) ;
%     title(UniqueCellTypes{uc})
% end
% 
% figure % by cell type histograms
% numColumns = 3 ;
% for uc = 1:length(UniqueCellTypes) ;
%     subplot(ceil(length(UniqueCellTypes)/numColumns),numColumns,uc)
%     plot(last_peak_HistBins,last_peak_byCellType_HistValues{uc}) ;
%     title(UniqueCellTypes{uc})
% end
% 
% figure % tuning curves for each cell type
% for uc = 1:length(UniqueCellTypes) ;
%     subplot(ceil(length(UniqueCellTypes)/numColumns),numColumns,uc)
%     for c=1:length(cell_i{uc}) ;
%         ci = cell_i{uc}(c) ;
%         errorbar(log10(Irel_unique),psth_peak_delta_trls_mean(:,ci),psth_peak_delta_trls_std(:,ci),psth_peak_delta_trls_std(:,ci))
%         hold on
%         title(UniqueCellTypes{uc})
%     end
% end
%     
% figure % average psth across all cells
% for ts=1:UniqueStimNum ; % for each trigger set 
%     subplot(ceil(UniqueStimNum/3),3,ts)
%     plot(psthTime,psth_mean(ts,:),'k','linewidth',3) ;
%     hold on
%     %plot(psthTime,psth{ds},'c') ;
%     plot(psthTime,psth_mean(ts,:)+psth_std(ts,:),'b:') ;
%     plot(psthTime,psth_mean(ts,:)-psth_std(ts,:),'b:') ;
%     %text(.9,.9,['NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),',fdur ',num2str(DataBlock(DB).DfParams.Ftime(ds))],'units', 'norm') 
% end
% 
% figure
% for c=1:NumCells ;
%     
%     subplot(1,2,1)
%     for ts=1:UniqueStimNum ; % for each trigger set 
%         %subplot(length(NDFs),length(Ftimes),plot_i(ds))
%         plot(psthTime,psth_delta{ts}(c,:)) ;
%         hold on
%         %xlim([-.2, .75])
%     end
%     hold off
%     
%     subplot(1,2,2)
%     errorbar(log10(Irel_unique),psth_peak_delta_trls_mean(:,c),psth_peak_delta_trls_std(:,c),psth_peak_delta_trls_std(:,c))
%     hold on 
%     if ~isnan(Threshold(c));
%         plot(Threshold(c),psth_peak_delta_trls_mean(Threshold(c)==log10(Irel_unique),c),'ro')
%     end
%     plot(log10(Irel_unique),psth_delta_peak(:,c),'g')
%     hold off
%     
%     title(num2str(c))
%     pause
% end
% 
% figure % example cell type individual cell data
% for c=1:length(cell_i{ExampleCellTypei}) ;
%     ci = cell_i{ExampleCellTypei}(c) ;
%     subplot(1,2,1)
%     for ts=1:UniqueStimNum ; % for each trigger set 
%         %subplot(length(NDFs),length(Ftimes),plot_i(ds))
%         plot(psthTime,psth_delta{ts}(ci,:)) ;
%         hold on
%         xlim([-.2, 1])
%     end
%     hold off
% 
%     subplot(1,2,2)
%     errorbar(log10(Irel_unique),psth_peak_delta_trls_mean(:,ci),psth_peak_delta_trls_std(:,ci),psth_peak_delta_trls_std(:,ci))
%     hold on 
%     if ~isnan(Threshold(ci));
%         plot(Threshold(ci),psth_peak_delta_trls_mean(Threshold(ci)==log10(Irel_unique),ci),'ro')
%     end
%     plot(log10(Irel_unique),psth_delta_peak(:,ci),'g')
%     hold off
% 
%     title(num2str(ci))
%     pause
% end 

% example cell data
if ~isfield(DataBlock(DB),'ExampleCelli') 
    ExampleCelli = input('Provide example cell # to export: ') ;
else
    if isempty(DataBlock(DB).ExampleCelli)
        ExampleCelli = input('Provide example cell # to export: ') ;
    else
        ExampleCelli = DataBlock(DB).ExampleCelli ;
    end
end
% 
% figure % look at rasters/psth for one cell for each flash
% for ts=1:UniqueStimNum ; % for each trigger set
%     subplot(10,1,1)
%     plot(psthTime,psth_delta{ts}(ExampleCelli,:)) ;
%     
%     subplot(10,1,2:10)
%     for t=1:length(full_trigger_set_i{ts}) ;
%         spk = dataRun.spikes{ExampleCelli}-dataRun.triggers(full_trigger_set_i{ts}(t)) ;
%         spk = spk(spk>=psthTime(1) & spk<=psthTime(end)) ;
%         plot(spk,ones(1,length(spk))*t,'k.')
%         axis([psthTime(1),psthTime(end),0,length(full_trigger_set_i{ts})]) 
%         hold on
%     end
%     title(num2str(ts))
%     hold off
%     pause
% end
%  
% % flash i for example cell data
% if ~isfield(DataBlock(DB),'ExampleFlashi') 
%     Flashi = input('Provide flash # for to export: ') ;
% else
%     if isempty(DataBlock(DB).ExampleFlashi)
%         Flashi = input('Provide flash # for to export: ') ;
%     else
%         Flashi = DataBlock(DB).ExampleFlashi ;
%     end
% end

            
% for igor
    
ForIgor = struct() ;
% 
ForIgor.ForMatlab(DB).ExampleCelli = ExampleCelli ;
ForIgor.ForMatlab(DB).ExampleCellTypei = ExampleCellTypei ;
%ForIgor.ForMatlab(DB).ExampleFlashi = Flashi ;
ForIgor.ForMatlab(DB).HistCellNumber = NumCells ;
% 
% VecName = ['logIrel','Db',num2str(DB)] ;
% ForIgor = setfield(ForIgor,VecName,log10(Irel_unique)) ; % relative intensity vector
% ForIgor.ForMatlab(DB).Irel_unique = Irel_unique ;
% 
% VecName = ['RspMean','Cell',num2str(ExampleCelli),'Db',num2str(DB)] ;
% ForIgor = setfield(ForIgor,VecName,psth_peak_delta_trls_mean(:,ExampleCelli)) ; % mean of trials
% ForIgor.ForMatlab(DB).psth_peak_delta_trls_mean = psth_peak_delta_trls_mean ;
% 
% VecName = ['RspStd','Cell',num2str(ExampleCelli),'Db',num2str(DB)] ;
% ForIgor = setfield(ForIgor,VecName,psth_peak_delta_trls_std(:,ExampleCelli)) ; % std of trials
% ForIgor.ForMatlab(DB).psth_peak_delta_trls_std = psth_peak_delta_trls_std ;
% 
% VecName = ['psthTime','Db',num2str(DB)] ;
% ForIgor = setfield(ForIgor,VecName,psthTime) ; % psth time vector
% ForIgor.ForMatlab(DB).psthTime = psthTime ;
% 
% for ts=1:UniqueStimNum ; % for each trigger set 
%     VecName = ['psth','Irel',num2str(ts),'Cell',num2str(ExampleCelli),'Db',num2str(DB)] ;
%     ForIgor = setfield(ForIgor, VecName, psth_delta{ts}(ExampleCelli,:)) ; % psth vector
%     ForIgor.ForMatlab(DB).psth_delta{ts} = psth_delta{ts} ;
% end
% 
% for t=1:length(full_trigger_set_i{Flashi}) ; % for each trial for 1 flash in 1 cell 
%     spk = dataRun.spikes{ExampleCelli}-dataRun.triggers(full_trigger_set_i{Flashi}(t)) ;
%     spk = spk(spk>=psthTime(1) & spk<=psthTime(end)) ;
%     spkY = ones(1,length(spk))*t ;
%     
%     VecName = ['Raster',num2str(t),'Irel',num2str(Flashi),'Cell',num2str(ExampleCelli),'Db',num2str(DB)] ;
%     ForIgor = setfield(ForIgor, VecName, spk) ; % spike times
%     
%     VecName = ['RasterY',num2str(t),'Irel',num2str(Flashi),'Cell',num2str(ExampleCelli),'Db',num2str(DB)] ;
%     ForIgor = setfield(ForIgor, VecName, spkY) ; % raster trial
%     
%     %ForIgor.ForMatlab(DB).Raster{t} = spk ;
% end

if ~isnan(ExampleCellTypei) ;% example cell type histograms
    VecName = ['ThreshHistValues', ExampleCellType,'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,Threshold_byCellType_HistValues{ExampleCellTypei}) ; % threshold hist values
    ForIgor.ForMatlab(DB).Threshold_byCellType_HistValues = Threshold_byCellType_HistValues ;
    
    VecName = ['FracNoThresh', ExampleCellType,'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,FractNoThresh_byCellType{ExampleCellTypei}) ; % threshold 
    ForMatlab(DB).Threshold_FracNoThresh_byCellType = FractNoThresh_byCellType ;

%     VecName = ['PkHistValues', ExampleCellType,'Db',num2str(DB)] ;
%     ForIgor = setfield(ForIgor,VecName,last_peak_byCellType_HistValues{ExampleCellTypei}) ; % hist values
%     ForIgor.ForMatlab(DB).last_peak_byCellType_HistValues = last_peak_byCellType_HistValues ;
% 
%     VecName = ['PkTimeHistValues', ExampleCellType,'Db',num2str(DB)] ;
%     ForIgor = setfield(ForIgor,VecName,last_peak_time_byCellType_HistValues{ExampleCellTypei}) ; % threshold hist values
%     ForIgor.ForMatlab(DB).last_peak_time_byCellType_HistValues = last_peak_time_byCellType_HistValues ;
% 
%     VecName = ['WidthHistValues', ExampleCellType,'Db',num2str(DB)] ;
%     ForIgor = setfield(ForIgor,VecName,last_halfWidth_byCellType_HistValues{ExampleCellTypei}) ; % threshold hist values
%     ForIgor.ForMatlab(DB).last_halfWidth_byCellType_HistValues = last_halfWidth_byCellType_HistValues ;
else
    ForIgor.ForMatlab(DB).Threshold_byCellType_HistValues = [] ;
    ForMatlab(DB).Threshold_FracNoThresh_byCellType = [] ;
%     ForIgor.ForMatlab(DB).last_peak_byCellType_HistValues = [] ;
%     ForIgor.ForMatlab(DB).last_peak_time_byCellType_HistValues = [] ;
%     ForIgor.ForMatlab(DB).last_halfWidth_byCellType_HistValues = [] ;
end
    

% all cell histograms
VecName = ['ThreshHistBins','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,Threshold_HistBins) ; % threshold bins
ForIgor.ForMatlab(DB).Threshold_HistBins = Threshold_HistBins ;

VecName = ['ThreshHistValues','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,Threshold_HistValues) ; % threshold hist values
ForIgor.ForMatlab(DB).Threshold_HistValues = Threshold_HistValues ;

VecName = ['BgRateHistBins','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,Bg_rate_HistBins) ; % bg rate bins
ForIgor.ForMatlab(DB).Bg_rate_HistBins = Bg_rate_HistBins ;

VecName = ['BgRateHistValues','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,Bg_rate_HistValues) ; % bg rate hist values
ForIgor.ForMatlab(DB).Bg_rate_HistValues = Bg_rate_HistValues ;

VecName = ['FracNoThresh','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,FracNoThresh) ; % fraction no threshold
ForIgor.ForMatlab(DB).Threshold_FracNoThresh = FracNoThresh ;

% VecName = ['PkHistBins','Db',num2str(DB)] ;
% ForIgor = setfield(ForIgor,VecName,last_peak_HistBins) ; % bins
% ForIgor.ForMatlab(DB).last_peak_HistBins = last_peak_HistBins ;
% 
% VecName = ['PkHistValues','Db',num2str(DB)] ;
% ForIgor = setfield(ForIgor,VecName,last_peak_HistValues) ; % hist values
% ForIgor.ForMatlab(DB).last_peak_HistValues = last_peak_HistValues ;
% 
% VecName = ['PkTimeHistBins','Db',num2str(DB)] ;
% ForIgor = setfield(ForIgor,VecName,last_peak_time_HistBins) ; % threshold bins
% ForIgor.ForMatlab(DB).last_peak_time_HistBins = last_peak_time_HistBins ;
% 
% VecName = ['PkTimeHistValues','Db',num2str(DB)] ;
% ForIgor = setfield(ForIgor,VecName,last_peak_time_HistValues) ; % threshold hist values
% ForIgor.ForMatlab(DB).last_peak_time_HistValues = last_peak_time_HistValues ;
% 
% VecName = ['WidthHistBins','Db',num2str(DB)] ;
% ForIgor = setfield(ForIgor,VecName,last_halfWidth_HistBins) ; % threshold bins
% ForIgor.ForMatlab(DB).last_halfWidth_HistBins = last_halfWidth_HistBins ;
% 
% VecName = ['WidthHistValues','Db',num2str(DB)] ;
% ForIgor = setfield(ForIgor,VecName,last_halfWidth_HistValues) ; % threshold hist values
% ForIgor.ForMatlab(DB).last_halfWidth_HistValues = last_halfWidth_HistValues ;
