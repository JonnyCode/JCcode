% this simululation will examine how a k-means cluster approach deals with
% a simple set of linear responses

% JC 6/19/17

% simulated linear filters
FilterTime=[.001:.001:.3] ;    % time vector

tpeak = ones(1,4)*.05 ;   % time of peak
peakRise = ones(1,4)*.010 ;   % rise of peak

ttrough = ones(1,4)*.11 ;   % time of trough
troughDecay = ones(1,4)*.06 ;    % decay of trough

peakAmp = [1,-1,1,-1] ;    % sign and amplitude of peak
troughAmp = [.2,-.2,1,-1] ;   % sign and amplitude of trough

for cl=1:length(tpeak) ;
    Filter(cl,:) = simFilter(FilterTime,tpeak(cl),peakRise(cl),peakAmp(cl),ttrough(cl),troughDecay(cl),troughAmp(cl)) ;
end

% stimulus
StimTime = [.001:.001:10] ;
stim = normrnd(0,10,1,length(StimTime)) ;
stim = lowPassFilter(stim, 1000, 50) ;

% lin response
for cl=1:length(tpeak) ;
    lr(cl,:) = conv(stim,[zeros(1,length(Filter)),Filter(cl,:)],'same') ;
end

% lin respons + noise
NoiseStd = 0 ; % number std (relative to signal) of noise
NumTrials = 100 ;
for cl=1:length(tpeak) ; % for each cell type
    for t=1:NumTrials ; % for each trial
        r{t}(cl,:) =lr(cl,length(FilterTime)+1:end) + normrnd(0,NoiseStd*std(lr(cl,:)),1,length(lr(cl,length(FilterTime)+1:end))) ;
    end
end

rmat = cell2mat(r) ;
rstim = repmat([1:length(lr(cl,length(FilterTime)+1:end))],1,NumTrials) ; % vector of response stimuli

% stim ensemble
for se = 1:length(StimTime)-length(FilterTime) ; 
    stimEnsem(se,:) = stim(se:se+length(FilterTime)) ;
end
%[stimEnsemPC,stimEnsemPrj,stimEnsemV] = pca(stimEnsem) ; % stimplify stim space - not used

% cluster stim and response
numClusters = [10,25] ;

for k=1:length(numClusters) ; % for each k in k-means
    stimClst{k} = kmeans(stimEnsem,numClusters(k)) ;

    resClst{k} = kmeans(rmat',numClusters(k)) ;
end

% joint distribution (column-->count of response cluster given a stim cluster)    
for k=1:length(numClusters) ; % for each k in k-means
    JointDist{k}= zeros(numClusters(k)) ; % initialize Joint Distribution
    for clstr = 1:numClusters(k) ; % for each cluster
        clstri = find(stimClst{k} == clstr) ; % time bins within stim cluster 
        for pnt = 1:length(clstri) ; % for each time bin in the stim cluster
            JointDist{k}(:,clstr) = JointDist{k}(:,clstr) + ...
                hist(resClst{k}(rstim==clstri(pnt)),[1:numClusters(k)])' ; % add the histogram 
        end
    end
end

% precent correct for optimal ML discrimitor
for k=1:length(numClusters) ; % for each k in k-means
    PercentCorrect(k) = sum(max(JointDist{k},[],2))/sum(JointDist{k}(:)) ;
end



