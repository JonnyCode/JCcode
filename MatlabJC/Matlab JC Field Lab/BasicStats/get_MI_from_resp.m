function I=get_MI_from_resp(spikes,bins,word_size,starts)

%inputs: spikes are organized by bins and trials
%        bins are trial_start:bin_size:trial_stop, for eg 0:.005:10
%        word_size is in number of bins, generally between 1 and 4
%        starts is a vector of start times for each trial

%outputs: I is information in bits


%p(wordi|t) for each t
maxnsp=max(max(spikes))+1; %maximum number of spikes in a bin
num_poss_words=maxnsp^word_size; %number of possible words
T=length(bins)-word_size+1; %number of time windows to average over in sliding time window

%count up how many times each word happens
wordi_counts=nan(T,num_poss_words);
for t=1:T
    wordi_over_trials=spikes(:,t:(t+word_size-1));
    
    %get word id for each trial
    word_id=nan(length(starts),1);
    for trial=1:length(starts)
        id=0;
        for sb=1:word_size
            id=id+wordi_over_trials(trial,sb)*maxnsp^(sb-1);  %gives unique id for each word       
        end
        word_id(trial)=id;
    end
    wordi_counts(t,:)=histcounts(word_id,0:num_poss_words,'Normalization','probability');
end

%p(wordi): above averaged for all t
prob_wordi=mean(wordi_counts,1);

%signal entropy
Hsig_i=-prob_wordi.*log2(prob_wordi);
Hsig=nansum(Hsig_i);

%noise entropy
Hnoise_t=nan(1,T);
for t=1:T
    temp=-wordi_counts(t,:).*log2(wordi_counts(t,:)); 
    %use log2(0)=0?
    Hnoise_t(t)=nansum(temp);
end
Hnoise=mean(Hnoise_t);

I=Hsig-Hnoise;

end