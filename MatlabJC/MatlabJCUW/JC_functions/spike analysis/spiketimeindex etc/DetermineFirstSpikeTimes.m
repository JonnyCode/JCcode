function [MeanFirstSpikeProbability,AnalyzeBursts,FirstSpikeTime,FirstSpikeTimeHistograms,FirstSpikeTimeStandardDeviations,TempFirstSpikeTime,CandidateBurstStart]=DetermineFirstSpikeTimes(SpikeTrain,PSTH,sigma3)



% determine which bursts should be analyzed - the bursts which are preceded by
% 3sigma(ms) of no activity) 

stimpoints=length(PSTH);

CandidateBurst=zeros(1,stimpoints);

for L=(sigma3*10)+1:stimpoints;                                     % ignore the first points (because you don't want to look back for the presence of a spike at negative times)
    if PSTH(L)>=1;                                  % when there is a non-zero value in the PSTH....
        if PSTH((L-round(sigma3*10)):L-1)<1;           % ....if that value is preceded by sigma3 ms of inactivity 
           CandidateBurst(L)=PSTH(L);        % ....then keep that timepoint
        end
    end
end

CandidateBurstStart=find(CandidateBurst>=1);


% now determine the time of the first spike following the identified periods of
% inactivity - spike must occur within simga3 ms of CandidateBurstStart



[b,c]=size(SpikeTrain)

FirstSpikeProbability=zeros(b,stimpoints);
TempFirstSpikeTime=zeros(b,stimpoints);

for m=1:b;                                                          % for each of the epochs
    for n=1:length(CandidateBurstStart);                                 % for each of the candidate burst periods in an epoch             
        range=CandidateBurstStart(n):CandidateBurstStart(n)+(sigma3*10);
        [r,x,z]=find(SpikeTrain(m,range)==1,1);
        TempFirstSpikeTime(m,CandidateBurstStart(n)+(x-1))=1;
        clear r x z range;
    end
end



% now determine, for each identified burst of activity, whether there exists a 
% spike within 3sigma following the the period of inactivity in at least
% 90% of the trials. 


for y=1:b
    for z=1:length(CandidateBurstStart);
        count=1;
        for x=CandidateBurstStart(z):CandidateBurstStart(z)+(sigma3*10)
            if TempFirstSpikeTime(y,x)>=1;
            FirstSpikeProbability(y,CandidateBurstStart(z))=count;
            count=count+1;
            end
        end
    end
end

        
MeanFirstSpikeProbability=mean(FirstSpikeProbability,1);
AnalyzeBursts=find(MeanFirstSpikeProbability>=0.9);


% now, for each burst of activity that results in a spike with sigma3 ms of
% the burst onset in >=90% of the epochs, get the first spike time 


for s=1:b;
    for t=1:length(AnalyzeBursts);
        range=AnalyzeBursts(t):AnalyzeBursts(t)+(sigma3*10);
        [u,v,w]=find(SpikeTrain(s,range)==1,1);
        if isempty(v)==1
            FirstSpikeTime(s,t)=0
        else
            FirstSpikeTime(s,t)=AnalyzeBursts(t)+(v-1);
        end
            clear u v w range;
     end
end

for n=1:length(AnalyzeBursts);
    if any(FirstSpikeTime(:,n)>0)
        temp=FirstSpikeTime(:,n);
        FirstSpikeTimeHistograms{n}=temp(temp>0)/10;                            % convert back to ms
        FirstSpikeTimeStandardDeviations(n)=std(FirstSpikeTimeHistograms{n})    % spike times converted back to ms by nature of previous line
    end
end


