function ForIgor = SquareAnalysisConcatWithBwMapping(DataBlock, DB, Params)

% adapted form 'SquareAnalysisConcat.m'.  Added stuff to include BWN
% srf/trf.  Also, included several additional paramaterizations of 

% JC 8/22/2016

% parameters
numTriggerBins = 4 ; % number of bins between each trigger 
ctStdThresh = 3 ; % std above the mean 
oldStimCodeFlag=1 ; % if dg stim params need to be adjusted
saveFigFlag = false ;

% load individual data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1) ;

for st = 1:4  ; % assumes data00smallSquares-bigSquare-ffPulse-driftingGrating

    TempPath = [DataBlock(DB).DataConcat(1:end-21),num2str(st),'/data00',num2str(st)] ; % first dataset is small square path
    StimPath{st} = [TempPath(1:end-15),'stimuli/s',TempPath(end-1:end),'.mat'] ;
    
    dataRunTemp = load_data(TempPath,opt) ;
    duration(st) = dataRunTemp.duration ;
    
    clear dataRunTemp
end

% analyze concatinated data set
dataRun = load_data(DataBlock(DB).DataConcat,opt) ;
sampleRate = dataRun.sampling_rate ;
    
% num cells
numCells = length(dataRun.spikes) ;

% time for each dataset
timeBands = cumsum([0,duration]) ;

% analyze each data set within concat data
for st=[1:4] ; % for each data set

    % time of spikes and triggers for this data set
    timeBand = [timeBands(st)+1/dataRun.sampling_rate, timeBands(st+1)] ;

    if st==1 || st==2 ;% analyze square data set

        % load stimulus
        load(StimPath{st}) ;

        % stim parameters 
        Trigs = dataRun.triggers(dataRun.triggers>=timeBand(1) & dataRun.triggers<=timeBand(2)) ;
        numTrigs = length(Trigs) ; % number of triggers
        numSquares = length(stim_out.x_vertices_all) ; % number of unique squares

        % srf
        % psth time bins 
        triggerDiff = diff([Trigs',Trigs(end)+max(diff(Trigs))]) ; % trigger interval
        for t=1:numTrigs ;
            psthTimeBins((t-1)*numTriggerBins+1:t*numTriggerBins+1) = ...
                [Trigs(t):triggerDiff(t)/numTriggerBins:Trigs(t)+triggerDiff(t)] ;
        end

        % spike rate
        for cl=1:numCells ; % for each cell
            spikes{cl} = dataRun.spikes{cl}(dataRun.spikes{cl}>=timeBand(1) & dataRun.spikes{cl}<=timeBand(2)) ;
            psth(cl,:) = histc(spikes{cl},psthTimeBins) ; % spike counts between bins
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

        % estimate a srf contrours from stvMax 

        for cl=1:numCells ; % for each cell
            
            maxV = max(stvMax{cl}(:)) ; % greatest value in stvMax
            [maxR,maxC] = find(stvMax{cl}==maxV) ; % max point in stv

            stvThreshTemp = mean(stvMax{cl}(:)) + std(stvMax{cl}(:))*ctStdThresh ; % x std above mean
            stvThresh = min(maxV,stvThreshTemp) ; % take the max point if the thresh is too high
            ct = contourc(stvMax{cl},[1,1]*std(stvMax{cl}(:))*ctStdThresh) ; % contours

            ctStarti = find(ct(1,:) == ct(1,1) ); % each individual ct
            numCt = length(ctStarti); % number of contours
            if numCt>0 ;
                ctStarti(end+1) = length(ct)+1 ;
            end

            cnt = 1 ; % counter
            for pg=1:numCt ; % for each polygon
                tempPoly = ct(1:2,ctStarti(pg)+1:ctStarti(pg+1)-1) ; % band
                if sum(inpolygon(maxC,maxR,tempPoly(1,:),tempPoly(2,:)))>0 ; % if the max points are within this polygon
                    srf{cl}{cnt} = tempPoly ; % designate the polygon a srf
                    cnt=cnt+1 ;
                end
            end  
        end

        % spike response

        clear psth*

        stimDuration = dataRun.triggers(3)-dataRun.triggers(1) ;
        for cl=1:numCells ; % for each cell
            for s=1:numSquares ; % for each square

                tr=find(stim_out.trial_list==s) ; % the square trials
                triggs = Trigs(tr*2-1) ; % ON triggers for each square

                [psth(s,:),psthTime] = get_smooth_psth(spikes{cl},triggs','stop',stimDuration) ;
                psth_max(s) = max(psth(s,:)) ;
            end

            [temp,maxi] = max(psth_max) ;
            psth_maxtrial(cl,:) = psth(maxi,:) ;
            psth_mean(cl,:) = mean(psth,1) ; % mean across all squares
            
            % max trial psth parameters
            OnPnts = round(mean(triggerDiff)*sampleRate) ;
            psth_maxtrial_params{cl}.On.peakAmp = max(psth_maxtrial(cl,1:OnPnts)) ;
            psth_maxtrial_params{cl}.On.mean = mean(psth_maxtrial(cl,1:OnPnts)) ;
            psth_maxtrial_params{cl}.On.peakDivMean = psth_maxtrial_params{cl}.On.peakAmp/psth_maxtrial_params{cl}.On.mean ;
            psth_maxtrial_params{cl}.On.halfPeakTime = find(psth_maxtrial(cl,1:OnPnts)>...
                (psth_maxtrial_params{cl}.On.peakAmp+psth_maxtrial_params{cl}.On.mean)/2, 1, 'first')*1/sampleRate ;
            
            psth_maxtrial_params{cl}.Off.peakAmp = max(psth_maxtrial(cl,OnPnts+1:end)) ;
            psth_maxtrial_params{cl}.Off.mean = mean(psth_maxtrial(cl,OnPnts+1:end)) ;
            psth_maxtrial_params{cl}.Off.peakDivMean = psth_maxtrial_params{cl}.Off.peakAmp/psth_maxtrial_params{cl}.Off.mean ;
            psth_maxtrial_params{cl}.Off.halfPeakTime = find(psth_maxtrial(cl,OnPnts+1:end)>...
                (psth_maxtrial_params{cl}.Off.peakAmp+psth_maxtrial_params{cl}.Off.mean)/2, 1, 'first')*1/sampleRate ;
            
            psth_maxtrial_params{cl}.On.OnOffIndex = (psth_maxtrial_params{cl}.On.peakAmp - psth_maxtrial_params{cl}.Off.peakAmp)/...
                (psth_maxtrial_params{cl}.On.peakAmp + psth_maxtrial_params{cl}.Off.peakAmp) ;
                
            
            tr=find(stim_out.trial_list==maxi) ; % the square trials
            triggs = Trigs(tr*2-1) ; % ON triggers for each square
            for trg=1:length(triggs) ; % for each trial
                spikeTimes{cl}{trg} = spikes{cl}(spikes{cl}>=triggs(trg) & spikes{cl}<triggs(trg)+stimDuration)-triggs(trg) ; % spike times
            end
            clear psth
        end
        
        % make neuron struct
        for cl=1:numCells ; % for each cell
            if st==1 ;
                NeuronStruct(cl).smallSquare.StvOn = stvOn{cl} ;
                NeuronStruct(cl).smallSquare.StvOff = stvOff{cl} ;
                NeuronStruct(cl).smallSquare.StvMax = stvMax{cl} ;
                
                NeuronStruct(cl).smallSquare.StaOn = staOn{cl} ;
                NeuronStruct(cl).smallSquare.StaOff = staOff{cl} ;
                NeuronStruct(cl).smallSquare.StaMax = staMax{cl} ;
                
                NeuronStruct(cl).smallSquare.Contours = srf{cl} ;
                
                NeuronStruct(cl).smallSquare.Spikes = spikeTimes{cl} ;
                NeuronStruct(cl).smallSquare.psth = psth_maxtrial(cl,:) ;
                NeuronStruct(cl).smallSquare.psth_params = psth_maxtrial_params{cl} ;
            elseif st==2 ;
                NeuronStruct(cl).bigSquare.StvOn = stvOn{cl} ;
                NeuronStruct(cl).bigSquare.StvOff = stvOff{cl} ;
                NeuronStruct(cl).bigSquare.StvMax = stvMax{cl} ;
                
                NeuronStruct(cl).bigSquare.StaOn = staOn{cl} ;
                NeuronStruct(cl).bigSquare.StaOff = staOff{cl} ;
                NeuronStruct(cl).bigSquare.StaMax = staMax{cl} ;
                
                NeuronStruct(cl).bigSquare.Contours = srf{cl} ;
                
                NeuronStruct(cl).bigSquare.Spikes = spikeTimes{cl} ;
                NeuronStruct(cl).bigSquare.psth = psth_maxtrial(cl,:) ;
                NeuronStruct(cl).bigSquare.psth_params = psth_maxtrial_params{cl} ;
            end
        end
    end
    clear psth* spike* srf 
    
    if st==3 ; % if stim is ff pulse
        
        Trigs = dataRun.triggers(dataRun.triggers>=timeBand(1) & dataRun.triggers<=timeBand(2)) ;
        stimDuration = 2 ;
        
        for cl=1:numCells ; % for each cell
            cnt = 1; % counter
            for trg=1:6:length(Trigs) ; % for each trial     
                spikeTimes{cl}{cnt} = dataRun.spikes{cl}(dataRun.spikes{cl}>=Trigs(trg) & ...
                    dataRun.spikes{cl}<Trigs(trg)+stimDuration)-Trigs(trg) ; % spike times
                cnt=cnt+1 ;
            end
            psth(cl,:) = get_smooth_psth(dataRun.spikes{cl},Trigs(1:3:end),'stop',stimDuration) ;
            
            % psth parameters
            OnPnts = round((Trigs(2)-Trigs(1))*sampleRate) ;
            psth_params{cl}.On.peakAmp = max(psth(cl,1:OnPnts)) ;
            psth_params{cl}.On.mean = mean(psth(cl,1:OnPnts)) ;
            psth_params{cl}.On.peakDivMean = psth_params{cl}.On.peakAmp/psth_params{cl}.On.mean ;
            psth_params{cl}.On.halfPeakTime = find(psth(cl,1:OnPnts)>...
                (psth_params{cl}.On.peakAmp+psth_params{cl}.On.mean)/2, 1, 'first')*1/sampleRate ;
   
            psth_params{cl}.Off.peakAmp = max(psth(cl,OnPnts+1:end)) ;
            psth_params{cl}.Off.mean = mean(psth(cl,OnPnts+1:end)) ;
            psth_params{cl}.Off.peakDivMean = psth_params{cl}.Off.peakAmp/psth_params{cl}.Off.mean ;
            psth_params{cl}.Off.halfPeakTime = find(psth(cl,OnPnts+1:end)>...
                (psth_params{cl}.Off.peakAmp+psth_params{cl}.Off.mean)/2, 1, 'first')*1/sampleRate ;
            
            psth_params{cl}.On.OnOffIndex = (psth_params{cl}.On.peakAmp - psth_params{cl}.Off.peakAmp)/...
                (psth_params{cl}.On.peakAmp + psth_params{cl}.Off.peakAmp) ;

            NeuronStruct(cl).ffPulse.spikes = spikeTimes{cl} ;
            NeuronStruct(cl).ffPulse.psth = psth(cl,:) ;
            NeuronStruct(cl).ffPulse.psth_params = psth_params{cl} ;
        end   
    end
    clear psth* spike*  
    
    if st==4 ; % if stim is drifting gratings
        TrialTrigInterval=10 ;
        
        % find dsi
        dataRunTemp = dataRun ;
        dataRunTemp.names.stimulus_path = [StimPath{st}(1:end-4),'.txt'] ;
        dataRunTemp.triggers  = dataRun.triggers(dataRun.triggers>=timeBand(1) & dataRun.triggers<=timeBand(2)) ;
        
        dataRunTemp = load_stim(dataRunTemp,'user_defined_trigger_interval', TrialTrigInterval) ;
                
        if oldStimCodeFlag==1 ;
            dataRunTemp = convert_stim_MatlabDG(dataRunTemp) ;
        end
        
        totaldur = 0 ;
        cell_ids = dataRunTemp.cell_ids ;
        [NumSpikesCell, StimComb] = get_spikescellstim(dataRunTemp,cell_ids,totaldur) ;

        [mag dsindex magmax magave angle rho theta num U V spave] = dscellanalysis(NumSpikesCell, StimComb) ;
        
        % find osi
        StimComb_u1 = unique(StimComb(:,1)) ; % unique spatial
        StimComb_u2 = unique(StimComb(:,2)) ; % unique temporal
        StimComb_u3 = unique(StimComb(:,3)) ; % unique directions
        ut = length(StimComb_u1)*length(StimComb_u2)*length(StimComb_u3) ; % total number of unique stim
        
        for stmTp=1:length(StimComb_u2) ; % for each unique temporal stim (this assumes length(StimComb_u1)=1)
            stm = find(ismember(StimComb(:,2),StimComb_u2(stmTp))==1) ;
            for cl=1:numCells 
                [Osi_vl(cl,stmTp),Osi_ah(cl,stmTp),Po_ah(cl,stmTp)] = OsiFinder(StimComb(stm,3)*pi/180, NumSpikesCell(cl,stm)) ;
                GrRsp(cl,stmTp) = mean(NumSpikesCell(cl,stm)) ; % mean rsp across directions
            end
        end

        for cl=1:numCells ;
            [thetaSort,thetai] = sort(theta{1}(cl,:)) ;
            NeuronStruct(cl).DriftGrat.fast.RadAngle = thetaSort ;
            NeuronStruct(cl).DriftGrat.fast.SpikesNum = rho{1}(cl,thetai) ;
            NeuronStruct(cl).DriftGrat.fast.dsi = dsindex{1}(cl) ;
            NeuronStruct(cl).DriftGrat.fast.VectorSum = mag{1}(cl) ;
            NeuronStruct(cl).DriftGrat.fast.VectorAngle = angle{1}(cl) ;
            NeuronStruct(cl).DriftGrat.fast.Osi_vl = Osi_vl(cl,1) ;
            NeuronStruct(cl).DriftGrat.fast.Osi_ah = Osi_ah(cl,1) ;
            NeuronStruct(cl).DriftGrat.fast.OrientPo = Po_ah(cl,1) ;
            NeuronStruct(cl).DriftGrat.fast.meanSpikeNum = GrRsp(cl,1) ;
            
            [thetaSort,thetai] = sort(theta{2}(cl,:)) ;
            NeuronStruct(cl).DriftGrat.slow.RadAngle = thetaSort ;
            NeuronStruct(cl).DriftGrat.slow.SpikesNum = rho{2}(cl,thetai) ;
            NeuronStruct(cl).DriftGrat.slow.dsi = dsindex{2}(cl) ; 
            NeuronStruct(cl).DriftGrat.slow.VectorSum = mag{2}(cl) ;
            NeuronStruct(cl).DriftGrat.slow.VectorAngle = angle{2}(cl) ;
            NeuronStruct(cl).DriftGrat.slow.Osi_vl = Osi_vl(cl,2) ;
            NeuronStruct(cl).DriftGrat.slow.Osi_ah = Osi_ah(cl,2) ;
            NeuronStruct(cl).DriftGrat.slow.OrientPo = Po_ah(cl,2) ;
            NeuronStruct(cl).DriftGrat.slow.meanSpikeNum = GrRsp(cl,2) ;
        end
    end      
end
 
% pick the best small contour
for cl = 1:numCells ;
    
    [riBig,ciBig] = find(NeuronStruct(cl).bigSquare.StvMax==max(NeuronStruct(cl).bigSquare.StvMax(:))) ; % location of max point of big square
    riBig = mean(riBig) ; % center X
    ciBig = mean(ciBig) ; % center Y
    
    [ri,ci] = find(NeuronStruct(cl).smallSquare.StvMax==max(NeuronStruct(cl).smallSquare.StvMax(:))) ; % location of max point of big square
    SqDist = sum(([ri,ci] - repmat([riBig,ciBig],length(ri),1)).^2,2) ; %euclidean distance between big and small square max spot locations   
    [m,mi] = min(SqDist) ;
    
    for c=1:length(NeuronStruct(cl).smallSquare.Contours) ;
        tempPoly = NeuronStruct(cl).smallSquare.Contours{c} ; % band
        if sum(inpolygon(ci(mi),ri(mi),tempPoly(1,:),tempPoly(2,:)))>0 ;
            NeuronStruct(cl).smallSquare.LikelyContour = NeuronStruct(cl).smallSquare.Contours{c} ; 
        end
    end
end
   
% cell type name default
for cl = 1:numCells ;
    NeuronStruct(cl).cellType = 'unidentified' ;
end

% binary white data run
% load Master data
dataRunBw = load_data(DataBlock(DB).BwPath{1}) ;
dataRunBw = load_neurons(dataRunBw) ;
dataRunBw = load_ei(dataRunBw, 'all') ;
dataRunBw = load_params(dataRunBw,'cell_type_depth', 5) ;
dataRunBw = load_sta(dataRunBw) ; % only necessary to get trf_time

% load receptive fields and filter
marks_params.thresh = 4.5;
dataRunBw = get_sta_summaries(dataRunBw, 'all','marks_params', marks_params);
 
filt_params.radius = 0.75;
dataRunBw = get_rfs_filtered(dataRunBw, 'all', 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');

% save spatial and temporal rf so can delete dataRunMaster later
for a=1:length(dataRunBw.stas.filt_rfs) ;
    srf{a} = dataRunBw.stas.filt_rfs{a}*dataRunBw.stas.polarities{a} ; % give the correct pix polarity 
end
trf = dataRunBw.stas.time_courses ; 

dataRunBw.stimulus.monitor_refresh = 60.35 ; % the monitor_refresh rate is not accurate
trf_time = -[0:dataRunBw.stas.depth-1]*dataRunBw.stimulus.interval/dataRunBw.stimulus.monitor_refresh ; 

% map using electrical images
cell_list_map = map_ei(dataRun, dataRunBw) ;

% save NeuronStruct    
save(['/Users/jcafaro/Desktop/Temp/Matfiles/SquareAnalysis',num2str(DB)],'NeuronStruct','-v7.3')


% figures
CellFig = figure ; % figure of each cell

for cl = 1:numCells ;
    figure(CellFig); clf
    
    colormap(brewermap([],'RdBu'))
    caxis([0,1]) 
    
    subplot(4,6,1:2)
    imagesc(NeuronStruct(cl).smallSquare.StvMax)
    hold on
    plot(NeuronStruct(cl).smallSquare.LikelyContour(1,:),NeuronStruct(cl).smallSquare.LikelyContour(2,:),'k')
    axis([0 400 0 400])
    title('small squares stv')
    
    subplot(4,6,3:4)
    imagesc(NeuronStruct(cl).bigSquare.StvMax)
    hold on
    plot(NeuronStruct(cl).bigSquare.Contours{1}(1,:),NeuronStruct(cl).bigSquare.Contours{1}(2,:),'k')
    axis([0 400 0 400])
    title('big squares stv')

    subplot(4,6,7:8)
    rasterPlot(NeuronStruct(cl).smallSquare.Spikes,[0,2])
    xlim([0,2])
    title('small squares trial with max response')
    
    subplot(4,6,9:10)
    rasterPlot(NeuronStruct(cl).bigSquare.Spikes,[0,2])
    xlim([0,2])
    title('big squares trial with max response')
    
    subplot(4,6,13:14)
    plot([1:length(NeuronStruct(cl).smallSquare.psth)]*1/sampleRate,NeuronStruct(cl).smallSquare.psth)
    xlim([0,2])
    ylim([0, max([NeuronStruct(cl).smallSquare.psth,NeuronStruct(cl).bigSquare.psth,NeuronStruct(cl).ffPulse.psth])])
    
    subplot(4,6,15:16)
    plot([1:length(NeuronStruct(cl).smallSquare.psth)]*1/sampleRate,NeuronStruct(cl).bigSquare.psth)
    xlim([0,2])
    ylim([0, max([NeuronStruct(cl).smallSquare.psth,NeuronStruct(cl).bigSquare.psth,NeuronStruct(cl).ffPulse.psth])])
    
    subplot(4,6,11:12)
    rasterPlot(NeuronStruct(cl).ffPulse.spikes,[0,2])
    xlim([0,2])
    title('full field step response')
    
    subplot(4,6,17:18)
    plot([1:length(NeuronStruct(cl).ffPulse.psth)]*1/sampleRate,NeuronStruct(cl).ffPulse.psth)
    xlim([0,2])
    ylim([0, max([NeuronStruct(cl).smallSquare.psth,NeuronStruct(cl).bigSquare.psth,NeuronStruct(cl).ffPulse.psth])])
    
    subplot(4,6,5)
    polar([NeuronStruct(cl).DriftGrat.slow.RadAngle,NeuronStruct(cl).DriftGrat.slow.RadAngle(1)],...
        [NeuronStruct(cl).DriftGrat.slow.SpikesNum,NeuronStruct(cl).DriftGrat.slow.SpikesNum(1)])
    hold on
    polar([NeuronStruct(cl).DriftGrat.slow.VectorAngle,NeuronStruct(cl).DriftGrat.slow.VectorAngle],...
        [0,NeuronStruct(cl).DriftGrat.slow.VectorSum],'r')
    title('slow drifting grating')
    
    subplot(4,6,6)
    polar([NeuronStruct(cl).DriftGrat.fast.RadAngle,NeuronStruct(cl).DriftGrat.fast.RadAngle(1)],...
        [NeuronStruct(cl).DriftGrat.fast.SpikesNum,NeuronStruct(cl).DriftGrat.fast.SpikesNum(1)])
    hold on
    polar([NeuronStruct(cl).DriftGrat.fast.VectorAngle,NeuronStruct(cl).DriftGrat.fast.VectorAngle],...
        [0,NeuronStruct(cl).DriftGrat.fast.VectorSum],'r')
    title('fast drifting grating')
    
    if ~isempty(cell_list_map{cl}) ; % if there is a mapped Bw cell
        Bwi = get_cell_indices(dataRunBw,cell_list_map{cl}) ; % get cell index
        
        subplot(4,6,19:20) % srf
        temp_rf = srf{Bwi} ;
        norm_rf = norm_image(temp_rf);
        imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
        set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])

        subplot(4,6,21:22) % temporal receptive field
        plot(fliplr(trf_time),trf{Bwi}) ;
    end
    
    if saveFigFlag ;
         print(gcf, '-djpeg', [saveFigPath,'Cell',num2str(cl)])
    else
        disp(['cell ',num2str(cl),': ', NeuronStruct(cl).cellType]) ;
        InputTxt = input('cell type:','s') ;
        if ~isempty(InputTxt)
            NeuronStruct(cl).cellType = InputTxt ;
        end
    end
end
      


CellFig = figure ;
for cl = 1:numCells ;
    figure(CellFig); clf
    
    subplot(4,6,1:2)
    imagesc(NeuronStruct(cl).smallSquare.StvMax)
    hold on
    plot(NeuronStruct(cl).smallSquare.LikelyContour(1,:),NeuronStruct(cl).smallSquare.LikelyContour(2,:),'k')
    axis([0 400 0 400])
    title('small squares stv')
    
    subplot(4,6,3:4)
    imagesc(NeuronStruct(cl).bigSquare.StvMax)
    hold on
    plot(NeuronStruct(cl).bigSquare.Contours{1}(1,:),NeuronStruct(cl).bigSquare.Contours{1}(2,:),'k')
    axis([0 400 0 400])
    title('big squares stv')

    subplot(4,6,7:8)
    rasterPlot(NeuronStruct(cl).smallSquare.Spikes,[0,2])
    xlim([0,2])
    title('small squares trial with max response')
    
    subplot(4,6,9:10)
    rasterPlot(NeuronStruct(cl).bigSquare.Spikes,[0,2])
    xlim([0,2])
    title('big squares trial with max response')
    
    subplot(4,6,13:14)
    plot([1:length(NeuronStruct(cl).smallSquare.psth)]*1/sampleRate,NeuronStruct(cl).smallSquare.psth)
    xlim([0,2])
    ylim([0, max([NeuronStruct(cl).smallSquare.psth,NeuronStruct(cl).bigSquare.psth,NeuronStruct(cl).ffPulse.psth])])
    
    subplot(4,6,15:16)
    plot([1:length(NeuronStruct(cl).smallSquare.psth)]*1/sampleRate,NeuronStruct(cl).bigSquare.psth)
    xlim([0,2])
    ylim([0, max([NeuronStruct(cl).smallSquare.psth,NeuronStruct(cl).bigSquare.psth,NeuronStruct(cl).ffPulse.psth])])
    
    subplot(4,6,11:12)
    rasterPlot(NeuronStruct(cl).ffPulse.spikes,[0,2])
    xlim([0,2])
    title('full field step response')
    
    subplot(4,6,17:18)
    plot([1:length(NeuronStruct(cl).ffPulse.psth)]*1/sampleRate,NeuronStruct(cl).ffPulse.psth)
    xlim([0,2])
    ylim([0, max([NeuronStruct(cl).smallSquare.psth,NeuronStruct(cl).bigSquare.psth,NeuronStruct(cl).ffPulse.psth])])
    
    subplot(4,6,5)
    polar([NeuronStruct(cl).DriftGrat.slow.RadAngle,NeuronStruct(cl).DriftGrat.slow.RadAngle(1)],...
        [NeuronStruct(cl).DriftGrat.slow.SpikesNum,NeuronStruct(cl).DriftGrat.slow.SpikesNum(1)])
    title('slow drifting grating')
    
    subplot(4,6,6)
    polar([NeuronStruct(cl).DriftGrat.fast.RadAngle,NeuronStruct(cl).DriftGrat.fast.RadAngle(1)],...
        [NeuronStruct(cl).DriftGrat.fast.SpikesNum,NeuronStruct(cl).DriftGrat.fast.SpikesNum(1)])
    title('fast drifting grating')
    
    if ~isempty(cell_list_map{cl}) ; % if there is a mapped Bw cell
        Bwi = get_cell_indices(dataRunBw,cell_list_map{cl}) ; % get cell index
        
        subplot(4,6,19:20) % srf
        temp_rf = srf{Bwi} ;
        norm_rf = norm_image(temp_rf);
        imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
        colormap(brewermap([],'RdBu'))
        caxis([0,1]) 
        set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])

        subplot(4,6,21:22) % temporal receptive field
        plot(fliplr(trf_time),trf{Bwi}) ;
    end
    
    disp(['cell ',num2str(cl),': ', NeuronStruct(cl).cellType]) ;
    InputTxt = input('cell type:','s') ;
    if ~isempty(InputTxt)
        NeuronStruct(cl).cellType = InputTxt ;
    end
end
            

CellTypesList = {} ;
cnt=1 ;
for cl = 1:numCells ;
    if ~ismember(NeuronStruct(cl).cellType,CellTypesList) ; % if this cell type is not in the list
        CellTypesList{cnt} = NeuronStruct(cl).cellType ;
        cnt = cnt+1 ;
    end
end
    
for typ=1:length(CellTypesList) ; % for each cell type
    figure
    for cl = 1:numCells ; % for each cell
        if strcmp(CellTypesList{typ},NeuronStruct(cl).cellType) ; % if it is this cell type
            plot(NeuronStruct(cl).smallSquare.LikelyContour(1,:),NeuronStruct(cl).smallSquare.LikelyContour(2,:)) ; % plot is small contour
            hold on
        end
    end
    title(CellTypesList{typ})
    axis([0 400 0 400])
end
    
    
    