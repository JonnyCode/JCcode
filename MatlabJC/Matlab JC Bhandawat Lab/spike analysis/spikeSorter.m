function [spikePntGroup] = spikeSorter(data, sampleRate, NegDiffThreshStd, numClusters, numPCs, plotFigs) 

% this function will sort spikes from scincillium recording into likely
% separate ORN groups.
% data - row vector
% sampleRate - Hz
% NegDiffThresh - delta mV
% numClusters - number of distinct groups to sort spikes into
% numPCs - number of principle components to sort by
% plotFigs - "true" ("false") if you want (don't want) figures ploted

% J. Cafaro 1/25/12

prePeakTime = .001 ;
postPeakTime = .004 ;
smthTime = .002 ;

prePeakPnts = prePeakTime*sampleRate ;
postPeakPnts = postPeakTime*sampleRate ;

% detect spikes
spikePnt = spikeFinderMuliUnit(data, sampleRate, NegDiffThreshStd, plotFigs) ;

% isolate spike waves
smthPnts = smthTime*sampleRate ;
dataHpf = highPassFilter(data,sampleRate,1) ; % high pass filter at 1 hz
dataHpfSmth = smooth(dataHpf,smthPnts) ;

upperBound = length(data)-postPeakPnts ;

spikePnt = spikePnt(spikePnt<=upperBound & spikePnt>prePeakPnts) ;
for a=1:length(spikePnt) ;
    spikeWave(a,:) = dataHpfSmth(spikePnt(a)-prePeakPnts:spikePnt(a)+postPeakPnts) ;
    spikeWave(a,:) = spikeWave(a,:)-mean(spikeWave(a,:)) ; % subtract mean
end

% run pca
for a=1:length(spikePnt) ; 
    spikeWaveDev(a,:) = spikeWave(a,:) - mean(spikeWave,1) ; % subtract sta
end
spikeWaveDevCov = cov(spikeWaveDev) ; % covariance matrix
[eigVecs,eigVals] = eig(spikeWaveDevCov) ; % eig values and vectors

% project major eigs onto spikewaves
for a=1:numPCs ;
    projValues(:,a) = spikeWave*eigVecs(:,end+1-a) ; 
end

% clustering
clusterInds = kmeans(projValues,numClusters) ; % cluster indicies

% pick out clustered spikes
for a=1:numClusters ;
    spikePntGroupTemp{a} = spikePnt(clusterInds==a) ;
    spikeGroupMeanTemp(a,:) = mean(spikeWave(clusterInds==a,:),1) ;
    spikeGroupStdTemp(a,:) = std(spikeWave(clusterInds==a,:),1);
end

% order clustered spikes by size
[temp,Maxi] = sort(max(spikeGroupMeanTemp,[],2)) ;
for a=1:numClusters ;
    spikePntGroup{a} = spikePntGroupTemp{Maxi(a)} ;
    spikeGroupMean(a,:) = spikeGroupMeanTemp(Maxi(a),:) ;
    spikeGroupStd(a,:) = spikeGroupStdTemp(Maxi(a),:) ;
end

% figures
if plotFigs ;
    colorVec = ['b','r','g','y','c'] ;

    figure
    plot(dataHpfSmth)
    hold on
    for a=1:numClusters ;
        plot(spikePntGroup{a},dataHpfSmth(spikePntGroup{a}),[colorVec(a),'*'])
    end
    xlabel('time')
    ylabel('data amplitude')

    figure
    plot(projValues(:,1),projValues(:,2),'*')
    hold on
    for a=1:numClusters ;
        plot(projValues(clusterInds==a,1),projValues(clusterInds==a,2),[colorVec(a),'o'])
    end
    xlabel('pc1')
    ylabel('pc2')

    figure
    for a=1:numClusters ;
        subplot(1,numClusters+1,a)
        plot([1:prePeakPnts+postPeakPnts+1],spikeWave(clusterInds==a,:)') 
        xlabel('time')
        ylabel('data amplitude')

        subplot(1,numClusters+1,numClusters+1)
        plot([1:prePeakPnts+postPeakPnts+1],spikeWave(clusterInds==a,:)',colorVec(a))
        hold on
        xlabel('time')
        ylabel('data amplitude') 
    end

    figure
    for a=1:numClusters ;
        plot([1:prePeakPnts+postPeakPnts+1],spikeGroupMean(a,:),colorVec(a))
        hold on
        plot([1:prePeakPnts+postPeakPnts+1],spikeGroupMean(a,:)+spikeGroupStd(a,:),[colorVec(a),':'])
        plot([1:prePeakPnts+postPeakPnts+1],spikeGroupMean(a,:)-spikeGroupStd(a,:),[colorVec(a),':'])
    end
    xlabel('time')
    ylabel('data amplitude')

    figure
    subplot(1,2,1)
    plot(diag(eigVals)/sum(diag(eigVals)),'*')
    xlabel('PC number')
    ylabel('fraction variance')

    subplot(1,2,2)
    plot([1:prePeakPnts+postPeakPnts+1],eigVecs(:,end-numPCs+1:end))
    xlabel('time')
    ylabel('eig vector amplitude')
end
