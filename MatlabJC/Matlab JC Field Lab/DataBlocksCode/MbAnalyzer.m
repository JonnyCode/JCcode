function ForIgor = MbAnalyzer(DataBlock, DB, Params)

% this function will analyze moving bar data

% JC 10/26/15
%DB=3; % TEMP


for MbType = 1:length(DataBlock(DB).MbPath) ; % for each moving bar group
    
    PhaseSpikeCountRsquared_HistX = [0:.1:1] ; 
    phase_width_HistX = [0:6:360] ; % degrees
    phase_dsi_HistX = [0:.1:1] ;
    phase_angle_HistX = [0:6:359] ; 
    widthRatio_HistX = [0:.02:2] ;
    dsiRatio_HistX = [0:.1:10] ; 

    
    % classify DS vs non-DS cells (using moving bars)
    opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1) ;
    dataRun = load_data(DataBlock(DB).MbPath{MbType},opt) ;
    dataRun.names.stimulus_path = [DataBlock(DB).MbPath{MbType}(1:end-15),'stimuli/s',DataBlock(DB).MbPath{MbType}(end-1:end),'.mat'] ;
    dataRun = load_stim_matlab(dataRun) ;
    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb_photons(dataRun,dataRun.cell_ids, dataRun.triggers(end)+6, 1) ;
    ds_struct = mbcellanalysis(NumSpikesCell, StimComb) ;

    a = 1; b = 2; % which parameters to use for classification
    Taylor1 = ds_struct.MAG{a,1}./((sum(ds_struct.RHO{a,1},2))');
    Taylor2 = ds_struct.MAG{b,1}./((sum(ds_struct.RHO{b,1},2))');
    Taylor_mean = (Taylor1+Taylor2)/2 ; % mean Taylor metric

    figure(1); clf ; 
    plot(log10(Taylor1), log10(Taylor2), '*')
    title('ndf0')
    xlabel('TP 1')
    ylabel('TP 2')
    hold on
    
    try
        load([DataBlock(DB).DSmatFilePath,'Db',num2str(DB),'Mbtype',num2str(MbType),'dsSelection.mat'])
        plot(log10(Taylor1(I)), log10(Taylor2(I)), 'ro')
    end
    
    if strcmp('n',input('use this selection (y/n)','s'))
        [x, y] = ginput;
        plot(x, y);
        IN = inpolygon(log10(Taylor1), log10(Taylor2), x, y);
        [~, I] = find(IN == 1);
        plot(log10(Taylor1(I)), log10(Taylor2(I)), 'ro')
        save([DataBlock(DB).DSmatFilePath,'Db',num2str(DB),'Mbtype',num2str(MbType),'dsSelection.mat'],'I') ;
    end
    ds_id = dataRun.cell_ids(I);

    % extract rates for low speed trials
    low_speed = 2;
    slow_indices = find(StimComb(:,2) == low_speed);

    [direction_list, sorted_dir_indices] = sort(StimComb(slow_indices,3), 'ascend');
    % sort indices into a ascending angles 
    num_reps = dataRun.stimulus.repetitions ;
    stim_duration = 10 ;
    
    % load skip cell indicies
    try
        load([DataBlock(DB).DSmatFilePath,'Db',num2str(DB),'Mbtype',num2str(MbType),'SkipRgc.mat']) ; 
    catch
        SkipRgc = zeros(1,length(ds_id)) ;
    end
        
    tuning_struct = [];
    analyzed_cell_num = 1;
    for rgc = 1:length(ds_id);

        temp_index = get_cell_indices(dataRun, ds_id(rgc));

        figure(1); clf ; 
        try
            load([DataBlock(DB).DSmatFilePath,'Db',num2str(DB),'Mbtype',num2str(MbType),'CellId',num2str(ds_id(rgc)),'phaseSelection.mat']) ; 
        catch
            clear raster_cut_times
        end

        for dr = 1:length(direction_list)
            temp_i = slow_indices(sorted_dir_indices(dr)) ;
            temp_triggers = find(dataRun.stimulus.trial_list==dataRun.stimulus.trial_list(temp_i)) ;
            
            temp_epochs = get_raster(dataRun.spikes{temp_index}, dataRun.stimulus.triggers(temp_triggers), 'stop', stim_duration,'plot', false);   
            mb_rasters(dr).rasters = temp_epochs;
            mb_rasters(dr).direction = direction_list(dr);

            subplot(3,4,dr)
            plot_raster(temp_epochs, 0, stim_duration)
            hold on
            
            if exist('raster_cut_times','var') ;
                plot([raster_cut_times(dr),raster_cut_times(dr)],[0,length(temp_triggers)],'r-')
            end

            if SkipRgc(rgc) == 1;
                title('SKIP')
            end

        end
        
        UserInStr = input('use this phase selection (y/n, skip)','s') ;
        if strcmp('n',UserInStr)
            raster_cut_times = [];
            [raster_cut_times,y] = ginput;
            
            save([DataBlock(DB).DSmatFilePath,'Db',num2str(DB),'Mbtype',num2str(MbType),'CellId',num2str(ds_id(rgc)),'phaseSelection.mat'],'raster_cut_times') ;
        end    
        if strcmp('skip',UserInStr) ;
            SkipRgc(rgc) = 1 ; 
        end
            
        if SkipRgc(rgc) == 0 ;
            for dr = 1:length(direction_list)
                temp_raster = mb_rasters(dr).rasters;
                spike_counter_one = 0;
                spike_counter_two = 0;
                for trial = 1:size(temp_raster,1)
                    temp_count = length(find(temp_raster{trial} < raster_cut_times(dr)));
                    spike_counter_one = temp_count + spike_counter_one;
                    temp_count = length(find(temp_raster{trial} > raster_cut_times(dr)));
                    spike_counter_two = temp_count + spike_counter_two;       
                end
                first_half_counter(dr) = spike_counter_one;
                second_half_counter(dr) = spike_counter_two;
            end

            ds_cell(rgc).first_half_counter = first_half_counter;
            ds_cell(rgc).second_half_counter = second_half_counter;

            % plot tuning curves for each phase of response
            figure(1); clf;
            subplot(2,1,1)
            polar([direction_list' 360] * (pi/180), [second_half_counter second_half_counter(1)], 'r')
            hold on
            polar([direction_list' 360] * (pi/180), [first_half_counter first_half_counter(1)], 'k')
            hold off

           % plot tuning curve
            subplot(2,1,2)
            plot(direction_list, second_half_counter ./ max(second_half_counter), 'ro')
            hold on
            plot(direction_list, first_half_counter ./ max(first_half_counter), 'ko')
            hold off
            drawnow
            pause
                
            PV = PolarVectorAddition([direction_list,first_half_counter']) ;
            first_phase_dir(analyzed_cell_num) = PV(1) ;
            first_phase_dsi(analyzed_cell_num) = PV(2)/sum(first_half_counter) ;
            
            PV = PolarVectorAddition([direction_list,second_half_counter']) ;
            second_phase_dir(analyzed_cell_num) = PV(1) ;
            second_phase_dsi(analyzed_cell_num) = PV(2)/sum(second_half_counter) ;
            
            temp_r = corrcoef(first_half_counter',second_half_counter') ;
            PhaseSpikeCountRsquared(analyzed_cell_num) = temp_r(1,2)^2 ;

            temp_width_one = circ_std(direction_list * (pi/180), first_half_counter', 30*pi/180) *180/pi;
            temp_width_two = circ_std(direction_list * (pi/180), second_half_counter', 30*pi/180) *180/pi;

            first_phase_width(analyzed_cell_num) = temp_width_one;
            second_phase_width(analyzed_cell_num) = temp_width_two;
            
            % figure for illustation (assumes 12 directions)
            figure
            subplot_place = [15,10,4,3, 2,6,11,16, 22,23,24,20] ;
            subplot(5, 5, [7:9,12:14,17:19])
            polar([direction_list' 360] * (pi/180), [first_half_counter+second_half_counter first_half_counter(1)+second_half_counter(1)], 'k')
            hold on
            polar([direction_list' 360] * (pi/180), [first_half_counter first_half_counter(1)], 'k--')
            polar([direction_list' 360] * (pi/180), [second_half_counter second_half_counter(1)], 'k:')
            polar([first_phase_dir(analyzed_cell_num),first_phase_dir(analyzed_cell_num)]* (pi/180),[0,first_phase_dsi(analyzed_cell_num)*sum(first_half_counter)],'k--')
            polar([second_phase_dir(analyzed_cell_num),second_phase_dir(analyzed_cell_num)]* (pi/180),[0,second_phase_dsi(analyzed_cell_num)*sum(second_half_counter)],'k:')
            for dr = 1:length(direction_list) ;
                subplot(5,5,subplot_place(dr))
                plot_raster(mb_rasters(dr).rasters, 0, stim_duration)
            end

            analyzed_cell_num = analyzed_cell_num + 1;
        end
    end
    save([DataBlock(DB).DSmatFilePath,'Db',num2str(DB),'Mbtype',num2str(MbType),'SkipRgc.mat'],'SkipRgc') ; 

    % population analysis  
    FirstSecondPhaseRatio_width = second_phase_width./first_phase_width ;
    FirstSecondPhaseRatio_dsi = second_phase_dsi./first_phase_dsi ;
    
    % histograms
    PhaseSpikeCountRsquared_Hist = hist(PhaseSpikeCountRsquared,PhaseSpikeCountRsquared_HistX) ; % r squared
    first_phase_width_Hist = hist(first_phase_width,phase_width_HistX) ; % phase width
    second_phase_width_Hist = hist(second_phase_width,phase_width_HistX) ; % phase width
    first_phase_dsi_Hist = hist(first_phase_dsi,phase_dsi_HistX) ; % phase vector dsi
    second_phase_dsi_Hist = hist(second_phase_dsi,phase_dsi_HistX) ; % phase vector dsi
    first_phase_dir_Hist = hist(first_phase_dir,phase_angle_HistX) ; % phase vector dir
    second_phase_dir_Hist = hist(second_phase_dir,phase_angle_HistX) ; % phase vector dir
    FirstSecondPhaseRatio_width_Hist = hist(FirstSecondPhaseRatio_width, widthRatio_HistX) ;
    FirstSecondPhaseRatio_dsi_Hist = hist(FirstSecondPhaseRatio_dsi, dsiRatio_HistX) ;

    first_phase_width_cHist = cumsum(first_phase_width_Hist)/sum(first_phase_width_Hist) ; % phase width
    second_phase_width_cHist = cumsum(second_phase_width_Hist)/sum(second_phase_width_Hist) ; % phase width
    first_phase_dsi_cHist = cumsum(first_phase_dsi_Hist)/sum(first_phase_dsi_Hist) ; % phase vector dsi
    second_phase_dsi_cHist = cumsum(second_phase_dsi_Hist)/sum(second_phase_dsi_Hist) ; %
    
    % delta dir histogram
    first_delta_dir_dist_median = 360 ;
    for DirFund = [0:1:90] ; % for each possible fundamental direction 
        delta_dir = [] ;
        for c = 1:length(first_phase_dir) ; % for each cell
            delta_dir(c) = min([acuteAngle(first_phase_dir(c),DirFund),...
            acuteAngle(first_phase_dir(c),DirFund+90),...
            acuteAngle(first_phase_dir(c),DirFund+180),...
            acuteAngle(first_phase_dir(c),DirFund+270)]) ; % the closest angle in this set of dir
        end
        delta_dir_median = median(delta_dir) ;
        if delta_dir_median<first_delta_dir_dist_median; % if the median is less than the smallest yet
            first_delta_dir_dist = delta_dir ;
            first_delta_dir_dist_median = delta_dir_median ;
        end
    end
    first_delta_dir_dist_Hist = hist(first_delta_dir_dist, phase_angle_HistX) ;
            
    second_delta_dir_dist_median = 360 ;
    for DirFund = [0:1:90] ; % for each possible fundamental direction 
        delta_dir = [] ;
        for c = 1:length(second_phase_dir) ; % for each cell
            delta_dir(c) = min([acuteAngle(second_phase_dir(c),DirFund),...
            acuteAngle(second_phase_dir(c),DirFund+90),...
            acuteAngle(second_phase_dir(c),DirFund+180),...
            acuteAngle(second_phase_dir(c),DirFund+270)]) ; % the closest angle in this set of dir
        end
        delta_dir_median = median(delta_dir) ;
        if delta_dir_median<second_delta_dir_dist_median; % if the median is less than the smallest yet
            second_delta_dir_dist = delta_dir ;
            second_delta_dir_dist_median = delta_dir_median ;
        end
    end              
    second_delta_dir_dist_Hist = hist(second_delta_dir_dist, phase_angle_HistX) ;
    
    % get tuning curves for analyzed cells 
    cell_num = 1 ;
    for rgc = 1:length(ds_id);
        if SkipRgc(rgc) == 0 ;
            first_half_count(cell_num,:) = ds_cell(rgc).first_half_counter ;
            second_half_count(cell_num,:) = ds_cell(rgc).second_half_counter ;
            
            cell_num = cell_num + 1 ;
        end
    end
            
    % get min firing rate of first and second phase spike count during thier prefered
    % stim
    minSpikeCount = min([max(first_half_count,[],2),max(second_half_count,[],2)],[],2)' ;
    
    % get degrees off of unity for direction of 1st and 2nd phases
    DegOffUnity = sqrt(2*(first_phase_dir-second_phase_dir).^2)/2 ;
    
    % select direction groups and measure their variance 
    figure
    for rgc=1:size(first_half_count,1) ;
        polar([first_phase_dir(rgc),first_phase_dir(rgc)] * (pi/180), [0,first_phase_dsi(rgc)*sum(first_half_count(rgc,:))],'b')
        hold on
        polar([second_phase_dir(rgc),second_phase_dir(rgc)] * (pi/180), [0,second_phase_dsi(rgc)*sum(second_half_count(rgc,:))],'r')
    axis tight
    end
    for drSelect=1:4 ;
        [x,y] = ginput ;
        Atemp = atan2d(y,x) ;
        Atemp(Atemp<0) = 360 + Atemp(Atemp<0) ;
        if max(Atemp)-min(Atemp)>acuteAngle(Atemp(1),Atemp(2))+2 ; % if its on the border
            drGroup_firstPhase{drSelect} = find(first_phase_dir>max(Atemp) | first_phase_dir<min(Atemp)) ;
            drGroup_secondPhase{drSelect} = find(second_phase_dir>max(Atemp) | second_phase_dir<min(Atemp)) ;
        else
            drGroup_firstPhase{drSelect} = find((first_phase_dir<max(Atemp)& first_phase_dir>min(Atemp))==1) ;
            drGroup_secondPhase{drSelect} = find((second_phase_dir<max(Atemp)& second_phase_dir>min(Atemp))==1) ;
        end
        for rgc=drGroup_firstPhase{drSelect} ;
            polar([first_phase_dir(rgc),first_phase_dir(rgc)] * (pi/180), [0,first_phase_dsi(rgc)*sum(first_half_count(rgc,:))],'k')
        end
        for rgc=drGroup_secondPhase{drSelect} ;
            polar([second_phase_dir(rgc),second_phase_dir(rgc)] * (pi/180), [0,second_phase_dsi(rgc)*sum(second_half_count(rgc,:))],'k')
        end

        if length(drGroup_firstPhase{drSelect})>2 ; % if their are enough cells to estimate a std
            drGroup_first_phase_std(drSelect) = circ_std(first_phase_dir(drGroup_firstPhase{drSelect})'*(pi/180))*(180/pi) ;
        else
            drGroup_first_phase_std(drSelect) = nan ;
        end
        if length(drGroup_secondPhase{drSelect})>2 ; % if their are enough cells to estimate a std
            drGroup_second_phase_std(drSelect) = circ_std(second_phase_dir(drGroup_secondPhase{drSelect})'*(pi/180))*(180/pi) ;
        else
            drGroup_second_phase_std(drSelect) = nan ;
        end
    end

    % get direction difference between first and second phase direction
    cell_num = 1 ;
    for rgc = 1:length(ds_id);
        if SkipRgc(rgc) == 0 ;
            PhaseDirDiff(cell_num) = acuteAngle(first_phase_dir(cell_num),second_phase_dir(cell_num)) ;
            cell_num = cell_num+1 ;
        end
    end
    PhaseDirDiff_Hist = hist(PhaseDirDiff,phase_angle_HistX) ;

    % save mat file of all data vectors
    save([DataBlock(DB).DSmatFilePath,'Db',num2str(DB),'Mbtype',num2str(MbType),'MbAnalysisMat.mat']) ;

    figure
    subplot(3,3,1)
    plot(first_phase_width, second_phase_width, 'ko')
    hold on
    plot([0:100],[0:100], 'k')
    axis([0 100 0 100])
    axis square
    xlabel('first phase tuning width')
    ylabel('second phase tuning width')
    hold off

    subplot(3,3,2)
    plot(first_phase_dir, second_phase_dir, 'ko')
    hold on
    plot([0:360],[0:360], 'k')
    axis([0 360 0 360])
    axis square
    xlabel('first phase direction')
    ylabel('second phase direction')
    
    subplot(3,3,3)
    plot(phase_width_HistX,first_phase_width_Hist,'-')
    hold on
    plot(phase_width_HistX,second_phase_width_Hist,'--')
    xlabel('phase tuning width')
    ylabel('number of observations')
    
    subplot(3,3,4)
    plot(phase_dsi_HistX,first_phase_dsi_Hist,'-')
    hold on
    plot(phase_dsi_HistX,second_phase_dsi_Hist,'--')
    xlabel('spike dsi')
    ylabel('number of observations')

    subplot(3,3,5)
    plot(phase_angle_HistX,first_phase_dir_Hist,'-')
    hold on
    plot(phase_angle_HistX,second_phase_dir_Hist,'--')
    xlabel('vector angle')
    ylabel('number of observations')
    
    subplot(3,3,6)
    plot(PhaseSpikeCountRsquared_HistX, PhaseSpikeCountRsquared_Hist, '-')
    hold on
    axis square
    xlabel('ON vs OFF r squared')
    ylabel('number of observations')
    
    subplot(3,3,7)
    plot(first_phase_dsi, first_phase_dir,'b*')
    hold on
    plot(second_phase_dsi, second_phase_dir,'r*')
    
    subplot(3,3,8)
    plot(widthRatio_HistX, FirstSecondPhaseRatio_width_Hist)
    
    subplot(3,3,9)
    plot(dsiRatio_HistX, FirstSecondPhaseRatio_dsi_Hist)
    
    figure
    plot(minSpikeCount,DegOffUnity,'*')
    
    figure
    plot(phase_angle_HistX, first_delta_dir_dist_Hist)
    hold on
    plot(phase_angle_HistX, second_delta_dir_dist_Hist)
    
    figure
    for rgc=1:size(first_half_count,1) ;
        subplot(5,ceil(size(first_half_count,1)/5),rgc)
        if max(first_half_count(rgc,:))>max(second_half_count(rgc,:))
            polar([direction_list', direction_list(1)] * (pi/180),[first_half_count(rgc,:),first_half_count(rgc,1)],'b')
            hold on
            polar([direction_list', direction_list(1)] * (pi/180),[second_half_count(rgc,:),second_half_count(rgc,1)],'r')
        else
             polar([direction_list', direction_list(1)] * (pi/180),[second_half_count(rgc,:),second_half_count(rgc,1)],'r')
             hold on
             polar([direction_list', direction_list(1)] * (pi/180),[first_half_count(rgc,:),first_half_count(rgc,1)],'b')
        end
  
        polar([first_phase_dir(rgc),first_phase_dir(rgc)] * (pi/180), [0,first_phase_dsi(rgc)*sum(first_half_count(rgc,:))],'b')
        polar([second_phase_dir(rgc),second_phase_dir(rgc)] * (pi/180), [0,second_phase_dsi(rgc)*sum(second_half_count(rgc,:))],'r')
        axis tight
    end
    
    figure
    for rgc=1:size(first_half_count,1) ;
        polar([first_phase_dir(rgc),first_phase_dir(rgc)] * (pi/180), [0,first_phase_dsi(rgc)*sum(first_half_count(rgc,:))],'b')
        hold on
        polar([second_phase_dir(rgc),second_phase_dir(rgc)] * (pi/180), [0,second_phase_dsi(rgc)*sum(second_half_count(rgc,:))],'r')
    axis tight
    end
    
    % For Igor
    ForIgor = struct() ; 
    
    VecName = ['firstPtuningW','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,first_phase_width) ; %
    
    VecName = ['secondPtuningW','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,second_phase_width) ; %
    
    
    VecName = ['firstPdir','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,first_phase_dir) ; %
    
    VecName = ['secondPdir','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,second_phase_dir) ; %

    
    VecName = ['firstPdsi','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,first_phase_dsi) ; %
    
    VecName = ['secondPdsi','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,second_phase_dsi) ; %
    
  
    VecName = ['phaseWidthHistX','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,phase_width_HistX) ; %
    
    VecName = ['firstPwidthHist','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,first_phase_width_Hist) ; %
    
    VecName = ['secondPwidthHist','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,second_phase_width_Hist) ; %
    
    VecName = ['firstPwidthHistc','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,first_phase_width_cHist) ; %
    
    VecName = ['secondPwidthHistc','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,second_phase_width_cHist) ; %
    
    
    VecName = ['phaseDsiHistX','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,phase_dsi_HistX) ; %
    
    VecName = ['firstPdsiHist','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,first_phase_dsi_Hist) ; %
    
    VecName = ['secondPdsiHist','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,second_phase_dsi_Hist) ; %
    
    VecName = ['firstPdsiHistc','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,first_phase_dsi_cHist) ; %
    
    VecName = ['secondPdsiHistc','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,second_phase_dsi_cHist) ; %
    
    
    VecName = ['phaseDirHistX','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,phase_angle_HistX) ; %
    
    VecName = ['firstPdirHist','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,first_phase_dir_Hist) ; %
    
    VecName = ['secondPdirHist','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,second_phase_dir_Hist) ; %
    
    VecName = ['firstDeltaDirDistHist','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,first_delta_dir_dist_Hist) ; %

    VecName = ['secondDeltaDirDistHist','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,second_delta_dir_dist_Hist) ; %
    
    VecName = ['phaseRsHistX','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,PhaseSpikeCountRsquared_HistX) ; %
    
    VecName = ['phaseRsHist','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,PhaseSpikeCountRsquared_Hist) ; % 
    
    
    VecName = ['widthRatioHistX','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,widthRatio_HistX) ; %
    
    VecName = ['widthRatioHist','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,FirstSecondPhaseRatio_width_Hist) ; % 
    
    
    VecName = ['dsiRatioHistX','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,dsiRatio_HistX) ; %
    
    VecName = ['dsiRatioHist','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,FirstSecondPhaseRatio_dsi_Hist) ; % 
    
    
    VecName = ['minSpikeCount','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,minSpikeCount) ; % 
    
    VecName = ['DegOffUnity','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,DegOffUnity) ; % 
    
    
    VecName = ['drGroupFirstPhaseStd','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,drGroup_first_phase_std) ; % 
    
    VecName = ['drGroupSecondPhaseStd','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,drGroup_second_phase_std) ; % 
    
    
    VecName = ['PhaseDirDiff','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,PhaseDirDiff) ; % 
    
    VecName = ['PhaseDirDiffHist','Db',num2str(DB),'Mbtype',num2str(MbType)] ;
    ForIgor = setfield(ForIgor,VecName,PhaseDirDiff_Hist) ; % 
    
    clearvars -except DataBlock DB Params
end


    
    
    
    