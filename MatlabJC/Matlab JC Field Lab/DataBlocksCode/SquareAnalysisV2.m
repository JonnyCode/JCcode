function ForIgor = SquareAnalysisV2(DataBlock, DB, Params)

% This function will analyze data from a moving flashing square
% modified from 'SquareAnalysis' to include mapping of EI, selection of
% specific region of stv and UNFINISHED!

% JC 11/7/2016

% parameters
numTriggerBins = 4 ; % number of bins between each trigger 
positionPlotFlag =1 ; % plot psth at center position of square

% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1) ;
dataRun = load_data(DataBlock(DB).MovingSquare,opt) ;

% num cells
numCells = length(dataRun.spikes) ;

% load stimulus
dataRun.names.stimulus_path = [DataBlock(DB).MovingSquare(1:end-15),'stimuli/s',DataBlock(DB).MovingSquare(end-1:end),'.mat'] ;
load(dataRun.names.stimulus_path) ;

% stim parameters 
numTrigs = length(dataRun.triggers) ; % number of triggers
numSquares = length(stim_out.x_vertices_all) ; % number of unique squares

%% srf
% psth time bins 
triggerDiff = diff([dataRun.triggers',dataRun.triggers(end)+max(diff(dataRun.triggers))]) ; % trigger interval
for t=1:numTrigs ;
    psthTimeBins((t-1)*numTriggerBins+1:t*numTriggerBins+1) = ...
        [dataRun.triggers(t):triggerDiff(t)/numTriggerBins:dataRun.triggers(t)+triggerDiff(t)] ;
end
        
% spike rate
for cl=1:numCells ; % for each cell
    psth(cl,:) = histc(dataRun.spikes{cl},psthTimeBins) ; % spike counts between bins
end

% spike rate change as function of time
for cl=1:numCells ; % for each cell
    psth_delta(cl,:) = psth(cl,:)-mean(psth(cl,:)) ;
    psth_var(cl,:) = (psth_delta(cl,:)).^2 ;
end

% stv and sta
ZeroMat = zeros(stim_out.x_end - stim_out.x_start + stim_out.stixel_width,...
        stim_out.y_end - stim_out.y_start + stim_out.stixel_width) ;

stimSpace = ZeroMat ;    

for cl=1:numCells ; % for each cell
    stvOn{cl}=ZeroMat ;
    stvOff{cl}=ZeroMat ;
    
    staOn{cl}=ZeroMat ;
    staOff{cl}=ZeroMat ;
    
    RespVarSum = sum(psth_var(cl,:)) ;
    RespSum = sum(psth(cl,:)) ;
    
    for s=1:numSquares ; % for each square position
        tr = stim_out.trial_list(s) ; % the square trial
        
        x_vert = stim_out.x_vertices_all{tr}-stim_out.x_start+1 ; % X-position
        y_vert = stim_out.y_vertices_all{tr}-stim_out.y_start+1 ; % Y-position

        tempMat=ZeroMat ;
        for q=1:length(x_vert(1,:)') ; % for each square (1 or 4) 
            tempMat([x_vert(1,q)'+1:x_vert(2,q)'],[y_vert(1,q)'+1:y_vert(3,q)'])=1 ; % make a binary square (the +1 makes point exclusive) 
        end
        
        if cl==1 ;
            stimSpace = stimSpace+tempMat ; % stim area 
        end
        
        Oni = s*2-1 ;
        RespVarOn = sum(psth_var(cl,(Oni-1)*numTriggerBins+1:Oni*numTriggerBins)) ; % response variance during square ON 
        RespOn = sum(psth(cl,(Oni-1)*numTriggerBins+1:Oni*numTriggerBins)) ; %
        
        Offi = s*2 ;
        RespVarOff = sum(psth_var(cl,(Offi-1)*numTriggerBins+1:Offi*numTriggerBins)) ; % response variance during square ON 
        RespOff = sum(psth(cl,(Offi-1)*numTriggerBins+1:Offi*numTriggerBins)) ; %
        
        stvOn{cl}=stvOn{cl}+tempMat*RespVarOn/RespVarSum ;
        stvOff{cl}=stvOff{cl}+tempMat*RespVarOff/RespVarSum ; 
        stvMax{cl}=reshape(max(stvOn{cl}(:),stvOff{cl}(:)),size(stvOn{cl}));
        
        staOn{cl}=staOn{cl}+tempMat*RespOn/RespSum ;
        staOff{cl}=staOff{cl}+tempMat*RespOff/RespSum ;  
        staMax{cl}=reshape(max(staOn{cl}(:),staOff{cl}(:)),size(staOn{cl}));
    end
end

stimSpace_binary = stimSpace==max(max(stimSpace)) ;

% plot ei % NOT SURE OF COORDINATES FOR IMAGESC(STA) 
load(DataBlock(2).TformEiPath) ; % load Tform (array-->monitor)
dataRun.piece.T_monitor_to_array = fliptform(Tform) ; % flip to (monitor-->array)
plot_ei(dataRun, dataRun.cell_ids(1), 'coordinate', 'monitor')

% estimate a srf contrours from stvMax 

for cl=1:numCells ; % for each cell
    
    ct = contourc(stvMax{cl},[1,1]*std(stvMax{cl}(:))*2) ; % contours
 
    ctStarti = find(ct(1,:) == ct(1,1) ); % each individual ct
    numCt = length(ctStarti); % number of contours
    if numCt>0 ;
        ctStarti(end+1) = length(ct)+1 ;
    end
   
    maxi = find(stvMax{cl}(:)==max(stvMax{cl}(:))) ; % max point in stv
    [maxR,maxC] = ind2sub(size(stvMax{cl}),maxi) ;
    
    cnt = 1 ; % counter
    for pg=1:numCt ; % for each polygon
        tempPoly = ct(1:2,ctStarti(pg)+1:ctStarti(pg+1)-1) ; % band
        if sum(inpolygon(maxC,maxR,tempPoly(1,:),tempPoly(2,:)))>0 ; % if the max points are within this polygon
            srf{cl}{cnt} = tempPoly ;
            cnt=cnt+1 ;
        end
    end  
end
  
%% spike response

clear psth*

if stim_out.sub_region == 1 ;
    numRows = stim_out.field_height * stim_out.shiftFactor / 2 ;
    numColumns = stim_out.field_width * stim_out.shiftFactor / 2 ;
else
    numRows = stim_out.field_height * stim_out.shiftFactor ;
    numColumns = stim_out.field_width * stim_out.shiftFactor ;
end
    
% trial number to subplot position
for s=1:numSquares ;
    Xcntr(s) = (stim_out.x_vertices_all{s}(1,1)+stim_out.x_vertices_all{s}(2,1))/2 ; % X-center
    Ycntr(s) = (stim_out.y_vertices_all{s}(1,1)+stim_out.y_vertices_all{s}(2,1))/2 ; % Y-center
end
[temp,temp,xi] = unique(Xcntr) ; % row 
[temp,temp,yi] = unique(Ycntr) ; % column 

% plot spike rate as function of spot center location
positionPlot = figure ;
psthPlot = figure ;

stimDuration = dataRun.triggers(3)-dataRun.triggers(1) ;
for cl=1:numCells ; % for each cell
    for s=1:numSquares ; % for each square
        
        tr=find(stim_out.trial_list==s) ; % the square trials
        triggs = dataRun.triggers(tr*2-1) ; % ON triggers for each square

        [psth(s,:),psthTime] = get_smooth_psth(dataRun.spikes{cl},triggs','stop',stimDuration) ;
        psth_var(s) = var(psth(s,:)) ;
    
        if positionPlotFlag ==1 ;
            figure(positionPlot)

            subploti = (xi(s)-1)*numColumns+yi(s) ;
            subplot(numRows,numColumns,subploti)
            plot(psthTime,psth(s,:)) ;
            set(gca,'ylim',[0,50])
            set(gca,'visible','off')
        end
    end
    
    [temp,maxi] = max(psth_var) ;
    psth_maxVar(cl,:) = psth(maxi,:) ;
    psth_mean(cl,:) = mean(psth) ;
    
    figure(psthPlot)
    clf
    plot(psthTime,psth_maxVar(cl,:)) ;
    hold on
    plot(psthTime,psth_mean(cl,:),'k') ;

    pause
end

%% figures

f1 = figure ;

for cl=1:numCells ; % for each cell
    figure(f1); clf
    subplot(2,3,1)
    imagesc(stvOn{cl}.*stimSpace_binary)
    colormap(brewermap([],'RdBu'))
    hold on
    colorbar 
    title(['cell: ',num2str(cl)])

    subplot(2,3,2)
    imagesc(stvOff{cl}.*stimSpace_binary)
    colormap(brewermap([],'RdBu'))
    hold on
    colorbar
    
    subplot(2,3,3)
    imagesc(stvMax{cl}.*stimSpace_binary)
    colormap(brewermap([],'RdBu'))
    colorbar
    hold on
    if ~isempty(srf{cl})
        plot(srf{cl}{1}(1,:),srf{cl}{1}(2,:),'k')
    end

    subplot(2,3,4)
    imagesc(staOn{cl}.*stimSpace_binary)
    colormap(brewermap([],'RdBu'))
    colorbar 

    subplot(2,3,5)
    imagesc(staOff{cl}.*stimSpace_binary)
    colormap(brewermap([],'RdBu')) 
    colorbar
    
    subplot(2,3,6)
    imagesc(staMax{cl}.*stimSpace_binary)
    colormap(brewermap([],'RdBu')) 
    colorbar
    pause  
end
    


