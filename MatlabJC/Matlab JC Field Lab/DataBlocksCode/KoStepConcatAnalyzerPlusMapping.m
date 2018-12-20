function ForIgor = KoStepConcatAnalyzerPlusMapping(DataBlock, DB, Params)

% modified version of KoStepConcatAnalyzer.m to include mapping to BW

% JC 3/4/2016

% parameters
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 
bl = 20 ; % number of gwgb pulses averaged in each block

DrugT = DataBlock(DB).FfPulseConcatDrugTimes ;

% load concatiniated data
dataRun = load_data(DataBlock(DB).FfPulseConcat) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun, 'all') ; % load electrical images

% load first set of concatinated data (assumes only first set of data may have been interupted)  
dataRunTemp = load_data(DataBlock(DB).FfPulse{1}) ;
dataRunTemp = load_neurons(dataRunTemp) ;

% find appropriate set of triggers
FirstSetExcess = rem(length(dataRunTemp.triggers),4) ;
trigs = [dataRun.triggers(4:length(dataRunTemp.triggers)-FirstSetExcess)',dataRun.triggers(length(dataRunTemp.triggers)+5:end)'] ;

blg = 4*bl ; % block group (total number of triggers per block)
num_bl = floor(length(trigs)/(blg)) ; % number of blocks 

for cl=1:length(dataRun.spikes) ; % for each cell
    for b = 1:num_bl ; % for each block of trials
        b_time(b) = trigs(b*blg-blg+1) ; % time of first data block trigger
        [psthTemp, binsTemp] = get_psth(dataRun.spikes{cl}, trigs(b*blg-blg+1:4:b*blg),'stop',12) ;
        psth{cl}(b,:) = psthTemp ;
    end 
    psth_mean(cl,:) = mean(psth{cl}) ; % mean psth
end
psth_time = binsTemp ;
 
BrkPnts = floor([0:size(psth_mean,2)/4:size(psth_mean,2)]) ; % points of stim change

for cl=1:length(dataRun.spikes) ; % for each cell
    [mx,mi] = max(psth_mean(cl,:)) ; % max of average psth 
    si = find(BrkPnts<mi,1,'last') ; % find stimulus that caused the max response

    [psth_peak(cl,:),mi] = max(psth{cl}(:,BrkPnts(si)+1:BrkPnts(si+1)),[],2) ;
    psth_rsp_mean(cl,:) = mean(psth{cl}(:,BrkPnts(si)+1:BrkPnts(si+1)),2) ;
    psth_rsp_trans(cl,:) = psth_peak(cl,:)./psth_rsp_mean(cl,:) ; % response transience (peak/mean)
    hp = psth_peak(cl,:)/2 ; % half peak
    for b = 1:num_bl ; % for each block of trials
        psth_duty(cl,b) = sum(psth{cl}(b,BrkPnts(si)+1:BrkPnts(si+1))>hp(b))/BrkPnts(2) ; % fraction of points above half peak
        psth_peak_time(cl,b) = psth_time(mi(b)) ; % (s) time post stim change of peak  
    end
end

% averages across cells
psth_peak_mean = mean(psth_peak,1) ;
psth_rsp_mean_mean = mean(psth_rsp_mean,1) ;
psth_rsp_trans_mean = mean(psth_rsp_trans,1) ;
psth_duty_mean = mean(psth_duty,1) ;
psth_peak_time_mean = mean(psth_peak_time,1) ;

psth_peak_std = std(psth_peak,0,1) ;
psth_rsp_mean_std = std(psth_rsp_mean,0,1) ;
psth_rsp_trans_std = std(psth_rsp_trans,0,1) ;
psth_duty_std = std(psth_duty,0,1) ;
psth_peak_time_std = std(psth_peak_time,0,1) ;


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

