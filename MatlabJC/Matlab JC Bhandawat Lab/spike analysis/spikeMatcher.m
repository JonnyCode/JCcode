function [spikePnt,usedSpikeTemplate] = spikeMatcher(data,projValueThresh,SpikeTemplate)

%spike detection algorthm which looks for waveforms within "data" that
%look like the "SpikeTemplate".  The waveform with the highest projection 
%becomes the "usedTemplate" for the entire trial

% jc 1/5/12

% spike template
SpikeTemplate = SpikeTemplate - mean(SpikeTemplate) ; 
lst = length(SpikeTemplate) ; 

% prep for speed
projValues = nan(1,length(data)-lst) ;
dataBlocks = nan(length(data)-lst,lst) ;

% project spike template
for a=1:length(data)-lst ;
    dataBlocks(a,:) = data(a:a+lst-1)-mean(data(a:a+lst-1)) ;
    projValues(a) = dataBlocks(a,:)*SpikeTemplate' ;
end

% find best match of spike template
[maxv,maxi] = max(projValues) ;
usedSpikeTemplate = dataBlocks(maxi,:) ;

% project new spike template
for a=1:length(data)-lst ;
    projValues(a) = dataBlocks(a,:)*usedSpikeTemplate' ;
end

% indentify spike time peaks
spikePnti = 0 ;
a = 1 ;
projThresh = projValueThresh*(usedSpikeTemplate*usedSpikeTemplate') ;
while a<=length(projValues) ; 

    if projValues(a)>projThresh ;
        spikePnti = spikePnti + 1 ;
        
        % local peak of projection values
        if length(projValues)-a+1>=lst ;
            [maxv,maxi1] = max(projValues(a:a+lst)) ;
        else
            [maxv,maxi1] = max(projValues(a:end)) ;
        end
        
        % local peak of data
        if length(data)-a+maxi1+1>=lst ;
            [maxv,maxi2] = max(data(a+maxi1:a+maxi1+lst)) ;
        else
            [maxv,maxi2] = max(data(a+maxi1:end)) ;
        end
        
        % spike point
        spikePnt(spikePnti) = a + maxi1 + maxi2 - 1 ;
 
        % set a
        if length(projValues)-spikePnt(spikePnti)+1>=lst ;
            a = spikePnt(spikePnti)+lst ;
        else
            a = length(projValues)+1 ;
        end
    else 
        a=a+1 ;
    end
end




