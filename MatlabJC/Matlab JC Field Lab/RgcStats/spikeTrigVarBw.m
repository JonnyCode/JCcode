function stv = spikeTrigVarBw(dataRun,cell_i,movie) 

% this function will get the spike triggered variance (stv), i.e. the
%variance (from a non-spike mean) of each movie pixel preceding a spike

% dataRun
% cell_i - cell indicies

% JC 12/30/15


% movie and params
frameRate = dataRun.stimulus.monitor_refresh/dataRun.stimulus.interval ; % frames/sec
%frameRate = size(movie,4)/dataRun.duration ; % frames/sec CALC IN THIS WAY DOES NOT WORK

mvRaw(:,:,:) = movie(:,:,1,:) ; % get rid of rgb for achromatic stim
mvRaw_mean = mean(mvRaw(:)) ; % average stim value
mvRaw_std = std(mvRaw(:)) ; % std of stim
mv = mvRaw - mvRaw_mean ; % subtract off average
mv_var= var(mv,0,3) ; % unbiased variance at each stixel

for c=1:length(cell_i) ; % for each cell
    sl = length(dataRun.stas.time_courses{cell_i(c)}) ; % spike triggered length
    
    % spike triggered average
    sta{c} = zeros(size(mv,1),size(mv,2),sl) ; % prep spike triggered variance
    spkFrames = ceil(dataRun.spikes{c}*frameRate) ;
    f0 = spkFrames(spkFrames>sl) ;
    for s=1:length(f0) ; % for each spike
        sta{c} = mv(:,:,(f0(s)-sl+1):f0(s))+sta{c} ;
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
    
    % SHOULD ALSO SUBTRACT OUT STIM VARIANCE BUT NEED TO CALC CORRECT STIM
    % VARIANCE
end
        
        
    
    
    