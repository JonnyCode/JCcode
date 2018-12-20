function [spk_times,peak_points] = spk_finder(raw_trace);

% detect times of negative going current spikes from cell attached pyramidal
% recordings during seizure (for Spain Lab rotation)
% J Cafaro 6/8/07


% set search parameters
sample_rate = .05 ;%ms          % sample rate of raw trace      
spk_thresh = -70 ;%mV          % negative spike threshold
ref = 1/sample_rate ;%points    % absolute refractory period

% find all points below a threshold
spk_points = find(raw_trace < spk_thresh);

% find breaks between points points crossing threshold (possible different
% spikes)
bpoints = find(diff(spk_points)>ref)' ;   % find where points are not consecuative and longer than refractory period
endpoints = [bpoints,length(spk_points)] ;
startpoints = [1,bpoints+1] ;

% find peak time of spikes
for a = 1:length(startpoints)  %for each spike...
    [peak(a),peak_point(a)] = min(raw_trace(spk_points(startpoints(a)):spk_points(endpoints(a)))) ; % lowest point on spike
end
    peak_points = peak_point + (spk_points(startpoints)-ones(length(startpoints),1))'; % points
    spk_times = peak_points * sample_rate ; % time of peak in ms

    

    

