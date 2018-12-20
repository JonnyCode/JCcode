function ForIgor = BwRepeatDrugAnalyzer(DataBlock, DB, Params)

% modified 'KoBwRepeatConcatAnalyzerPlusMapping.m' for just bw repeat data

% JC 11/15/2018


% data set select
RunAsScript = true ;
if RunAsScript ;
    DB = 57 ; % 
    [DataBlock,Params] = DataBlocks_KO ;
end

BwRepeatPathNum = 1 ;

% parameters
Color_list = {'k','r','b','g'} ; % order of colors for each 

DrugT = DataBlock(DB).BwRepeatDrugTimes ;

% load data
dataRun = load_data(DataBlock(DB).BwRepeat{BwRepeatPathNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun, 'all') ; % load electrical images

NumCells = length(dataRun.spikes) ;

% ei centers
numElectrodeLayers = 2 ;
for cells = 1:NumCells ; % for each cell
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
for cells = 1:NumCells ; % for each cell
    for cnd =1:length(DrugTrials)-1 ; % for each condition 
        for t=1:(DrugTrials(cnd+1)-DrugTrials(cnd)) ; % for each trial
            [psthTemp,binsTemp]  = get_psth(dataRun.spikes{cells},AllTriggers(DrugTrials(cnd)+t-1),'stop',10) ;
            psth{cells}{cnd}(t,:) = psthTemp ;
        end
        psth_mean{cells}(cnd,:) = mean(psth{cells}{cnd},1) ;
        psth_std{cells}(cnd,:) = std(psth{cells}{cnd},[],1) ;
    end
end

% psth change relative to 1st condition
for cells = 1:NumCells ; % for each cell
    for cnd =2:length(DrugTrials)-1 ; % for each change in condition 
        psth_mean_dprime{cells}(cnd) = mean((psth_mean{cells}(1,:)-psth_mean{cells}(cnd,:))./ ...
            sqrt(((psth_std{cells}(1,:).^2+psth_std{cells}(cnd,:).^2)/2)+.001));
    end
end

%% cross correlation between psths (looking for duplicates)
CorrThreshold = .9 ;
 
psth_AllMean = nan(NumCells,length(binsTemp)) ; % prep mat
for cells = 1:NumCells ; % for each cell
    psth_AllMean(cells,:) = mean(psth_mean{cells},1) ;
end
    
psth_corrMat = corr(psth_AllMean') ; % correlation between all psths
 
% list all pairs that exceed CorrThreshold
[c1,c2] = find(psth_corrMat.*(eye(NumCells)==0)>CorrThreshold) ;




%% map location of cells on epiflourescence image
Tform = map_Ei_to_camera(DataBlock(DB).IrImagePath, DataBlock(DB).BwRepeat{BwRepeatPathNum}) ; % get transform to camera

%figure
% CheckOn = [46,47,64,65,82,92,108,110,161,162,178,182,212,228,255,258,263,287,289,...
%     294,322] ; % Db56
CheckOn = [86,111,147,162,184,186,240] ;

im_array = imread(DataBlock(DB).EpiImagePath);

figure
imshow(im_array);
hold on
cnd = 2 ;
for cells=1:NumCells ;
    if ismember(cells,CheckOn)
        col = 'r' ;
    else
        col = 'c' ;
    end
    
    ctrImage = tformfwd(Tform,EiCtr(cells,1),EiCtr(cells,2)) ;
    plot(ctrImage(1),ctrImage(2),'+','Color',col)
    %plot(ctrImage(1),ctrImage(2),'o','MarkerSize',100*abs(psth_mean_dprime{cells}(cnd))+.001,'Color',col)
end


%% figures

% plot psth and raster - faster

for cells = 1:length(dataRun.spikes) ; % for each cell
    figure(1)
    clf
    subplot(2,1,1)
    for a =1:length(DrugTrials)-1 ; % for each condition 
        spikes_by_trials = get_raster(dataRun.spikes{cells}, AllTriggers(DrugTrials(a):DrugTrials(a+1)-1),...
            'stop', 10,'tic_color',Color_list{a},'foa',-1, 'first_tic',DrugTrials(a),'dots_flag',true) ;
        hold on
    end
    axis([0 10 1 DrugTrials(end)])
    drawnow
    
    subplot(2,1,2)
    for a =1:length(DrugTrials)-1 ; % for each condition 
        plot(binsTemp,psth_mean{cells}(a,:),Color_list{a})
        hold on
        plot(binsTemp,psth_mean{cells}(a,:)+psth_std{cells}(a,:),[':',Color_list{a}])
        plot(binsTemp,psth_mean{cells}(a,:)-psth_std{cells}(a,:),[':',Color_list{a}])
    end
    title(num2str(cells))
    
    pause
    
end



% plot psth and raster - fast 
for cells = 1:length(dataRun.spikes) ; % for each cell
    figure(1)
    clf
    subplot(2,1,1)
    RasterImagePlotter(dataRun.spikes{cells}, AllTriggers,'LineTrials',DrugTrials)
    axis([0 10 1 DrugTrials(end)])
    drawnow
    
    subplot(2,1,2)
    for a =1:length(DrugTrials)-1 ; % for each condition 
        plot(binsTemp,psth_mean{cells}(a,:),Color_list{a})
        hold on
        plot(binsTemp,psth_mean{cells}(a,:)+psth_std{cells}(a,:),[':',Color_list{a}])
        plot(binsTemp,psth_mean{cells}(a,:)-psth_std{cells}(a,:),[':',Color_list{a}])
    end
    title(num2str(cells))
    pause    
end


% plot psth and raster - slow

for cells = 1:length(dataRun.spikes) ; % for each cell
    figure(1)
    clf
    subplot(2,1,1)
    for a =1:length(DrugTrials)-1 ; % for each condition 
        spikes_by_trials = get_raster(dataRun.spikes{cells}, AllTriggers(DrugTrials(a):DrugTrials(a+1)-1),...
            'stop', 10,'tic_color',Color_list{a},'foa',-1, 'first_tic',DrugTrials(a)) ;
    end
    axis([0 10 1 DrugTrials(end)])
    drawnow
    
    subplot(2,1,2)
    for a =1:length(DrugTrials)-1 ; % for each condition 
        plot(binsTemp,psth_mean{cells}(a,:),Color_list{a})
        hold on
        plot(binsTemp,psth_mean{cells}(a,:)+psth_std{cells}(a,:),[':',Color_list{a}])
        plot(binsTemp,psth_mean{cells}(a,:)-psth_std{cells}(a,:),[':',Color_list{a}])
    end
    title(num2str(cells))
    pause    
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



check = [] ;

