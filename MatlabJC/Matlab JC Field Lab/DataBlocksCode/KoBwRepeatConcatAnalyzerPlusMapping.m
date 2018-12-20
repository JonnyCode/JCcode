function ForIgor = KoBwRepeatConcatAnalyzerPlusMapping(DataBlock, DB, Params)

% modified version of KoStepConcatAnalyzerPlusMapping.m to use for bw
% repeat data instead of steps

% JC 7/25/2016

% parameters
Color_list = {'k','r','b','g'} ; % order of colors for each 

DrugT = DataBlock(DB).BwRepeatConcatDrugTimes ;

% load concatiniated data
dataRun = load_data(DataBlock(DB).BwRepeatConcat) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun, 'all') ; % load electrical images

% ei centers
numElectrodeLayers = 2 ;
for cells = 1:length(dataRun.spikes) ; % for each cell
    EiCtr(cells,:) = get_ei_com(dataRun, dataRun.cell_ids(cells), numElectrodeLayers) ;
end

% drug triggers (assumes 6 trials every 10 sec repeat)
AllTriggers = dataRun.triggers(1:6:end) ; % triggers at begining of each trial
DrugTimes = [0,DrugT,dataRun.duration] ; % initial times of condition changes

for a = 1:length(DrugTimes) ;
    [m, mi] = min(abs(AllTriggers - DrugTimes(a))) ; % nearest trial start time
    DrugTrials(a) = mi ;
end
    
% get psth
for cells = 1:length(dataRun.spikes) ; % for each cell
    for a =1:length(DrugTrials)-1 ; % for each condition 
        [psthTemp, binsTemp] = get_psth(dataRun.spikes{cells},AllTriggers(DrugTrials(a):DrugTrials(a+1)-1),'stop',10) ;
        psth{cells}(a,:) = psthTemp ;
    end
end

% psth change 1st condition
for cells = 1:length(dataRun.spikes) ; % for each cell
    for a =2:length(DrugTrials)-1 ; % for each condition
        psth_DivCntrl{cells}(a,:) = psth{cells}(a,:)./psth{cells}(1,:) ;
        psth_Mse(cells,a) = sum((psth{cells}(a,:)-psth{cells}(1,:)).^2) ;
    end
end

% Binary white data run
% load Master data
dataRunMaster = load_data(DataBlock(DB).BwPath{1}) ;
dataRunMaster = load_neurons(dataRunMaster) ;
dataRunMaster = load_ei(dataRunMaster, 'all') ;
dataRunMaster = load_params(dataRunMaster,'cell_type_depth', 5) ;

% get cell types 
for a = 1:length(dataRunMaster.cell_types) ;
    celltypes{a} = dataRunMaster.cell_types{a}.name ;
end

UniqueCellTypes = unique(celltypes) ;
if isempty(UniqueCellTypes{1}) ;
    UniqueCellTypes = UniqueCellTypes(2:end) ;
end
    
% map using electrical images
cell_list_map = map_ei(dataRunMaster, dataRun) ;

% cells ids in slave for each UniqueCellType set in master data
for uc = 1:length(UniqueCellTypes) ;
    Tempi = get_cell_indices(dataRunMaster, UniqueCellTypes{uc}) ; % cell indicy of master
    master_id{uc} = [] ;
    for a=1:length(Tempi) ;
        if ~isempty(cell_list_map{Tempi(a)}) ;
            master_id{uc} = [master_id{uc},dataRunMaster.cell_ids(Tempi(a))] ; % cell id of master
        end
    end
    cell_ids{uc} = cell2mat(cell_list_map(Tempi)) ; % id of slave
    if ~isempty(cell_ids{uc}) ;
        cell_i{uc} = get_cell_indices(dataRun, cell_ids{uc}) ; % indicy of slave
    end
end

% cell sta 
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

% NOT SET FOR BW DATA
% last blocks in each condition
b_cntl = find(DrugT(1)>b_time,1,'last') ; % control block
b_drug1 = find(DrugT(2)>b_time & 1400 >b_time,1,'last') ; % drug block
b_drug2 = find(DrugT(2)>b_time,1,'last') ; % drug block
b_wash = num_bl ; % wash block

% more parameters for last set data blocks organized by cell type
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    lc = length(cell_i{uc}) ; % number of cell of that type
    if lc>0 ;
        for cells=1:lc ; % for each cell of this type
            for bkp = 1:length(BrkPnts)-1 ; % for each stim (g,w,g,b)
                psth_cntrl_BkpPeak{uc}(cells,bkp) = max(psth{cell_i{uc}(cells)}(b_cntl,BrkPnts(bkp)+1:BrkPnts(bkp+1))) ;
                psth_cntrl_BkpMean{uc}(cells,bkp) = mean(psth{cell_i{uc}(cells)}(b_cntl,BrkPnts(bkp)+1:BrkPnts(bkp+1))) ;

                psth_drug1_BkpPeak{uc}(cells,bkp) = max(psth{cell_i{uc}(cells)}(b_drug1,BrkPnts(bkp)+1:BrkPnts(bkp+1))) ;
                psth_drug1_BkpMean{uc}(cells,bkp) = mean(psth{cell_i{uc}(cells)}(b_drug1,BrkPnts(bkp)+1:BrkPnts(bkp+1))) ;

                psth_drug2_BkpPeak{uc}(cells,bkp) = max(psth{cell_i{uc}(cells)}(b_drug2,BrkPnts(bkp)+1:BrkPnts(bkp+1))) ;
                psth_drug2_BkpMean{uc}(cells,bkp) = mean(psth{cell_i{uc}(cells)}(b_drug2,BrkPnts(bkp)+1:BrkPnts(bkp+1))) ;
                
                psth_wash_BkpPeak{uc}(cells,bkp) = max(psth{cell_i{uc}(cells)}(b_wash,BrkPnts(bkp)+1:BrkPnts(bkp+1))) ;
                psth_wash_BkpMean{uc}(cells,bkp) = mean(psth{cell_i{uc}(cells)}(b_wash,BrkPnts(bkp)+1:BrkPnts(bkp+1))) ;

                psth_bkp_BkpPeak_DeltaDrug1{uc}(cells,bkp) = DiffOverSum(psth_drug1_BkpPeak{uc}(cells,bkp),psth_cntrl_BkpPeak{uc}(cells,bkp)) ;
                psth_bkp_BkpMean_DeltaDrug1{uc}(cells,bkp) = DiffOverSum(psth_drug1_BkpMean{uc}(cells,bkp),psth_cntrl_BkpMean{uc}(cells,bkp)) ;
                
                psth_bkp_BkpPeak_DeltaDrug2{uc}(cells,bkp) = DiffOverSum(psth_drug2_BkpPeak{uc}(cells,bkp),psth_cntrl_BkpPeak{uc}(cells,bkp)) ;
                psth_bkp_BkpMean_DeltaDrug2{uc}(cells,bkp) = DiffOverSum(psth_drug2_BkpMean{uc}(cells,bkp),psth_cntrl_BkpMean{uc}(cells,bkp)) ;
            end
            % psth of organized by cell type and condition
            psth_cntrl{uc}(cells,:) = psth{cell_i{uc}(cells)}(b_cntl,:) ;
            psth_drug1{uc}(cells,:) = psth{cell_i{uc}(cells)}(b_drug1,:) ;
            psth_drug2{uc}(cells,:) = psth{cell_i{uc}(cells)}(b_drug2,:) ;
            psth_wash{uc}(cells,:) = psth{cell_i{uc}(cells)}(b_wash,:) ; 
        end

        psth_cntrl_mean(uc,:) = mean(psth_cntrl{uc},1) ;
        psth_drug1_mean(uc,:)= mean(psth_drug1{uc},1) ;
        psth_drug2_mean(uc,:)= mean(psth_drug2{uc},1) ;
        psth_wash_mean(uc,:) = mean(psth_wash{uc},1) ; 
    end
end

% figures

% plot psth and raster

for cells = 1:length(dataRun.spikes) ; % for each cell
    figure(1)
    clf
    for a =1:length(DrugTrials)-1 ; % for each condition 
        spikes_by_trials = get_raster(dataRun.spikes{cells}, AllTriggers(DrugTrials(a):DrugTrials(a+1)-1),...
            'stop', 10,'tic_color',Color_list{a},'foa',-1, 'first_tic',DrugTrials(a)) ;
    end
    axis([0 10 1 DrugTrials(end)])
    drawnow
    
    figure(2)
    clf
    for a =1:length(DrugTrials)-1 ; % for each condition 
        plot(binsTemp,psth{cells}(a,:),Color_list{a})
        hold on
    end
    title(num2str(cells))
    pause    
end

% psth change organized by cell type
ChangeFig = figure ;
for uc = 1:length(UniqueCellTypes) ;
    figure(ChangeFig)
    clf
    for a =1:length(DrugTrials)-1 ; % for each condition 
        for cells = 1:length(cell_i{uc}) ; % for each cell of this type
            plot(binsTemp,psth_DivCntrl{cell_i{uc}(cells)}(a,:),Color_list{a})
            hold on
        end
    end
    
    drawnow
    title([UniqueCellTypes{uc},': ',num2str(cell_i{uc})])
    pause    
    hold off
end
    
% psth change organized by cell type - change colapsed across time 
figure
for uc = 1:length(UniqueCellTypes) ;
    for a =2:length(DrugTrials)-1 ; % for each condition 
        for cells = 1:length(cell_i{uc}) ; % for each cell of this type
            plot(uc,psth_Mse(cell_i{uc}(cells),a)./mean(psth{cell_i{uc}(cells)}(a,:))^2,['*',Color_list{a}])
            hold on
        end
    end
end

% psth change organized by cell type and location
figure
for uc = 1:length(UniqueCellTypes) ;
    %figure
    for a =3; % for each condition 
        for cells = 1:length(cell_i{uc}) ; % for each cell of this type
            [X,Y] = drawEllipse([[EiCtr(cell_i{uc}(cells),1),EiCtr(cell_i{uc}(cells),2)] 10 10]) ;
            plot(X,Y,'k','lineWidth',10*psth_Mse(cell_i{uc}(cells),a)./max(psth_Mse(cell_i{uc}(:),a))) 
           %plot(EiCtr(cell_i{uc}(cells),1),EiCtr(cell_i{uc}(cells),2),'o','MarkerSize',psth_Mse(cell_i{uc}(cells),a)./mean(psth{cell_i{uc}(cells)}(a,:))^2)
            hold on
        end
    end
end

% mse (does magintude of relative drug change depend on array location or cell type?)
figure

subplot(2,1,1)
for c = 1:length(dataRunMaster_cell_ids) ; % for each with a BW mapped ref
    if ~isempty(cell_list_map{c}) ; % if there is a mapped cell
        lw = MseRelDrug_DivMseRelWashLfMean(slave_c(c))/max(MseRelDrug_DivMseRelWashLfMean) ; 
        
        if ~isnan(lw) ;
            [X,Y] = drawEllipse([ctr{c} rad{c} angle{c}]) ;
            if ~any(isnan([X,Y])) ;
                [X,Y] = tformfwd(coord_tform, X, Y) ;
                plot(X,Y,'k','color',[1-lw,1-lw,1-lw],'linewidth',lw*3)
                hold on
            end
        end
    end
end   


% BELOW NOT UPDATED TO BW
figure
for cl=1:length(dataRun.spikes) ; % for each cell
    clf
    
    subplot(4,1,1)
    plot(psth_time,psth_mean(cl,:),'c','linewidth',4)
    hold on
    %plot(psth_time,psth{cl})
    for b = 1:num_bl ; % for each block of trials
        if DrugT(1)<b_time(b) && b_time(b)<DrugT(2) ; 
            plot(psth_time,psth{cl}(b,:),'r','linewidth',2*(1.1-b/num_bl))
        elseif DrugT(1)>b_time(b) ;
            plot(psth_time,psth{cl}(b,:),'k','linewidth',2*(1.1-b/num_bl))
        elseif DrugT(2)<b_time(b)
            plot(psth_time,psth{cl}(b,:),'b','linewidth',2*(1.1-b/num_bl))
        end    
        hold on
    end
    ylabel('firing rate (hz)')
    xlabel('time (s)')

        
    subplot(4,1,2)
    [ax,h1,h2] = plotyy(b_time,psth_peak(cl,:),b_time,psth_rsp_mean(cl,:)) ;
    h1.Marker = '*' ; h1.LineStyle = 'none' ; h2.Marker = '*' ; h2.LineStyle = 'none' ;
    hold on
    plot(DrugT,[psth_peak(cl,1),psth_peak(cl,1)],'r')
    ylabel('peak rate (hz)')
    xlabel('block start time (s)')
    
    subplot(4,1,3)
    plot(b_time,psth_peak_time(cl,:),'*')
    hold on
    plot(DrugT,[psth_peak_time(cl,1),psth_peak_time(cl,1)],'r')
    ylabel('peak time (s)')
    xlabel('block start time (s)')
    
    subplot(4,1,4)
    [ax,h1,h2] = plotyy(b_time,psth_duty(cl,:),b_time,1./psth_rsp_trans(cl,:)) ;
    h1.Marker = '*' ; h1.LineStyle = 'none' ; h2.Marker = '*' ; h2.LineStyle = 'none' ;
    hold on
    plot(DrugT,[psth_duty(cl,1),psth_duty(cl,1)],'r')
    ylabel('duty (%)')
    xlabel('block start time (s)')
    
    pause
end   


figure
subplot(3,1,1)
[ax,h1,h2] = plotyy(b_time,psth_peak_mean,b_time,psth_rsp_mean_mean) ;
h1.Marker = '*' ; h1.LineStyle = 'none' ; h2.Marker = '*' ; h2.LineStyle = 'none' ;
hold on
plot(DrugT,[psth_peak_mean(1),psth_peak_mean(1)],'r')
ylabel('peak rate (hz)')
xlabel('block start time (s)')

subplot(3,1,2)
plot(b_time,psth_peak_time_mean,'*')
hold on
plot(DrugT,[psth_peak_time_mean(1),psth_peak_time_mean(1)],'r')
ylabel('peak time (s)')
xlabel('block start time (s)')

subplot(3,1,3)
[ax,h1,h2] = plotyy(b_time,psth_duty_mean,b_time,1./psth_rsp_trans_mean) ;
h1.Marker = '*' ; h1.LineStyle = 'none' ; h2.Marker = '*' ; h2.LineStyle = 'none' ;
hold on
plot(DrugT,[psth_duty_mean(1),psth_duty_mean(1)],'r')
ylabel('duty (%)')
xlabel('block start time (s)')

% compare step responses for each cell type
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
 
        plot(psth_time,psth{cell_i{uc}(cells)}(b_cntl,:),'k')
        hold on
        plot(psth_time,psth{cell_i{uc}(cells)}(b_drug1,:),'r:')
        plot(psth_time,psth{cell_i{uc}(cells)}(b_drug2,:),'r')
        plot(psth_time,psth{cell_i{uc}(cells)}(b_wash,:),'b')
        
        title(['slave i:',num2str(cell_i{uc}(cells)),' master id:',num2str(master_id{uc}(cells))])
    end
end

% all cell together with parameters
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    if lc>0 ;
        for cells=1:lc ; % for each cell of this type
            subplot(5,1,4)% delta peak
            plot([1:length(BrkPnts)-1],psth_bkp_BkpPeak_DeltaDrug{uc}(cells,:),'color',[.8,.8,.8]) ;
            hold on

            subplot(5,1,5)% delta peak
            plot([1:length(BrkPnts)-1],psth_bkp_BkpMean_DeltaDrug{uc}(cells,:),'color',[.8,.8,.8]) ;
            hold on
        end
    %     subplot(5,1,1:3)% single example
    %     plot(psth_time,psth{cell_i{uc}(1)}(b_cntl,:),'color',[.8,.8,.8])
    %     hold on
    %     plot(psth_time,psth{cell_i{uc}(1)}(b_drug,:),'color',[.8,.6,.6])
    %     plot(psth_time,psth{cell_i{uc}(1)}(b_wash,:),'color',[.6,.6,.8])

        % means
        subplot(5,1,1:3)% 
        plot(psth_time, psth_cntrl_mean(uc,:),'k')
        hold on
        plot(psth_time, psth_drug1_mean(uc,:),'r:')
        plot(psth_time, psth_drug2_mean(uc,:),'r')
        plot(psth_time, psth_wash_mean(uc,:),'b')
        xlabel('time (s)')
        ylabel('spike rate (hz)')

        subplot(5,1,4)% delta peak
        plot([1:length(BrkPnts)-1],mean(psth_bkp_BkpPeak_DeltaDrug{uc},1),'k*-','linewidth',2) ;
        plot([1:length(BrkPnts)-1],zeros(1,length(BrkPnts)-1),'r:')
        xlabel('stim (g,w,g,b')
        ylabel('Change Peak ')

        subplot(5,1,5)% delta peak
        plot([1:length(BrkPnts)-1],mean(psth_bkp_BkpMean_DeltaDrug{uc},1),'k*-','linewidth',2) ;
        plot([1:length(BrkPnts)-1],zeros(1,length(BrkPnts)-1),'r:')
        xlabel('stim (g,w,g,b')
        ylabel('Change Mean ')
    end
end


check = [] ;

