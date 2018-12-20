function [number_exacts, number_matches, percent_match_ms, percent_match_quartms]=nearest_neighbor(spk_rastor4_zero_data1,spk_rastor4_zero_data2) ;

% This function will check to assess what precentage of spikes from a pair
% are within a given time from each other.
% J Cafaro 6/19/07

sumrastor = spk_rastor4_zero_data1 + spk_rastor4_zero_data2 ; % add together spike rastors

i = find(sumrastor==2) ; % find those times when spikes happen at the same time
number_exacts = length(i); % the number of spikes that happen at the exact same time
sumrastor(i) = 1 ; % get rid of all twos

for a=0:19;  % for each distance were interested within 1ms (20points at 20kHz)....
dist=[1,zeros(1,a),1] ;   % number of points between spikes
matches{a+1} = find(conv(sumrastor,dist)==2); % indices of all events with two spikes separated by 2 x dist ms
number_matches(a+1) = length(matches{a+1}) ; % number of spikes that have neihboor within a points of each other
end

number_spks_data1 = sum(spk_rastor4_zero_data1) ;   %number of spikes in each rastor
number_spks_data2 = sum(spk_rastor4_zero_data2) ;
min_number_spks = min(number_spks_data1,number_spks_data2) ; % number spikes in trace w/ fewer spikes

total_matches_ms = sum(number_matches)+ number_exacts ;    % total number of spikes that have another spike within 1ms
percent_match_ms = total_matches_ms*100/min_number_spks ; % percentage of spikes that have another spike within 1ms

total_matches_quartms = sum(number_matches(1:5))+ number_exacts ;    % total number of spikes that have another spike within 1ms
percent_match_quartms = total_matches_quartms*100/min_number_spks ; % percentage of spikes that have another spike within 1ms

