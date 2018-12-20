function stv = spikeTrigVarBlur(dataRun,cell_i,movie, stixBlurr, frameBlurr) 

%this function will average over space and time before calculating the spike 
% triggered variance (stv)

% JC 1/5/15

% movie and params
frameRate = dataRun.stimulus.monitor_refresh/dataRun.stimulus.interval ; % frames/sec
FrameShift = floor(frameBlurr/2) ;
spatialBlurr = ones(frameBlurr) ;

mvRaw(:,:,:) = movie(:,:,1,:) ; % get rid of rgb for achromatic stim
mvRaw_mean = mean(mvRaw(:)) ; % average stim value
mv = mvRaw - mvRaw_mean ; % subtract off average

for c=1:length(cell_i) ; % for each cell
    sl = length(dataRun.stas.time_courses{cell_i(c)}) ; % spike triggered length
    spkFrames = ceil(dataRun.spikes{c}*frameRate) ; % spike frames
    f0 = spkFrames(spkFrames>sl+FrameShift & spkFrames<(size(mv,3)-FrameShift)) ; % subset of spike frames
    
    % blurr (sliding average over space and time)
    stvTemp = zeros(size(mv,1),size(mv,2),sl) ; % prep stv
     for s=1:length(f0) ; % for each spike
         
        % spatial average 
        for preFrame = 1:sl+2*FrameShift ;
            blurrFrame(:,:,preFrame) = conv2(mv(:,:,f0(s)+FrameShift-preFrame+1),spatialBlurr,'same') ;
        end

        % temporal average
        for preFrame = 1:sl ;
           blurrFrame2(:,:,preFrame) = mean(blurrFrame(:,:,preFrame:preFrame+frameBlurr-1),3) ;
        end
        
        stvTemp = stvTemp+blurrFrame2.^2 ;
        
     end
     stv{c} = stvTemp/length(f0) ;
end
     
        
        
        
    
    
    
    
    
    % spike triggered average
    sta{c} = zeros(size(mv,1),size(mv,2),sl) ; % prep spike triggered variance
    spkFrames = ceil(dataRun.spikes{c}*frameRate) ;
    f0 = spkFrames(spkFrames>30) ;
    for s=1:length(f0) ; % for each spike
        sta{c} = mv(:,:,f0(s)-29:f0(s))+sta{c} ;
    end
    sta{c} = sta{c}/length(f0) ; % sta
    
    % find sig stixels and make time course
    sta_std = std(sta{c}(:)) ;
    sta_mean = mean(sta{c}(:)) ;
    sta_stdAbove = abs(sta{c}-sta_mean)/sta_std ;
    sta_stdAbove_stix = max(sta_stdAbove,[],3) ;
    sta_wieghts = repmat(sta_stdAbove_stix,[1,1,sl]) ;
    tcTemp = mean(mean(sta{c}.*sta_wieghts.^4,2),1) ; % time course is weighted average of stixel significance raised to 4th power
    tc(1,:) = tcTemp(1,1,:) ;
    
    tcm = nans(size(sta{c})) ; % prep time course matrix
    for tcf = 1:sl ; % for each time course frame
       tcm(:,:,tcf) = tc(tcf) ; % populate matrix 
    end
    
    % find spike triggered variance
    stvTemp = zeros(size(mv,1),size(mv,2)) ; % prep spike triggered variance
    for s=1:length(f0) ; % for each spike
        af = mean(mv(:,:,f0(s)-sl+1:f0(s)).*tcm,3) ; % weighted average of frames
        stvTemp = stvTemp + af.^2 ; % square and add 
    end
    stv{c} = stvTemp/length(f0); % stv
    
    % blurred and squarred
     for s=1:length(f0) ; % for each spike
        for preFrame = 0:length(dataRun.stas.time_courses{c})-1 ;
            blurrFrame(:,:,1,preFrame+1) = conv2(mv(:,:,1,f0(s)-preFrame),spatialBlurr,'same') ;
        end
        
        for preFrame = 0:length(dataRun.stas.time_courses{c})-1 ;
        %    blurrFrame(:,:,1,preFrame+1) = mean(mv(:,:,1,f0(s)-preFrame-XXX) ;
        end
        
        sta{c} = mv(:,:,1,f0(s)-29:f0(s))+sta{c} ;
        af = mean(mv(:,:,1,f0(s)-29:f0(s)).*tc,4) ;
        stv{c} = stv{c} + af.^2 ; % square and add 
    end
    
    stv{c} = stv{c}/length(spks) ; % variance 
end
        
        
    
    
    