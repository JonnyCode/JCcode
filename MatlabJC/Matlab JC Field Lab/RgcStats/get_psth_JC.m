function psth = get_smooth_psth(spike_times, trigger_times, varargin)

% JC 4/24/15
% function will make a psth from a set of trigger times.  Modified from
% "get_psth" to make sliding time bins and not histogram.  Also fixed
% problems with "start" param.
%
% optional parameters
%   start       0           start will be trigger_times + start
%   stop        []          stop time will be trigger_times - stop
%   first_trial 1           allows user to start with a later trial
%   sample_rate 20000       sample rate of spike
%   plot_hist   false       whether to generate PSTH on the fly
%   bin_size    0.1        bins size of PSTH
%   hist_color  [0 0 0]         color of PSTH
%   foa         0          calls to set_up_figure_or_axes function,
%

% SET UP OPTIONAL ARGUMENTS
p = inputParser;
p.addParamValue('start', 0, @isnumeric);
p.addParamValue('stop', [], @isnumeric);
p.addParamValue('first_trial', 1, @isnumeric)
p.addParamValue('plot_hist', false, @islogical);
p.addParamValue('bin_size', 0.1, @isnumeric);
p.addParamValue('hist_color', [0 0 0]);
p.addParamValue('foa', 0);
p.addParamValue('labels', true, @islogical);

% parse inputs
p.parse(varargin{:});
params = p.Results;
    
% figure out the duration of each trial
if isempty(params.stop)
    params.stop = mean(diff(trigger_times));
end

% compute the number of trials to use
num_trials = length(trigger_times) - params.first_trial + 1;

% spikeTrain prep
psthTime = [start:1/sample_rate:stop] ; % time vector
SpikeTrain = zeros(length(trigger_times),length(psthTime)) ;

for trg = 1:length(trigger_times) ; % for each trigger
    spike_times_tc = spike_times - trigger_times(trg) ; % trigger centered
    SpikeTrain(trg,:) = 1*(psthTime == spike_times_tc) ; 
end
    
    
    


% parset spike times into binned trials to generate a psth
for i = params.first_trial:length(trial_begin_times)
    % shift spike times by trial frames
    ith_shifted_times = spike_times - trial_begin_times(i);
    % identify spike times within each trial
    ith_trial_times = ith_shifted_times(ith_shifted_times >= params.start & ith_shifted_times <= params.stop);
    % bin the spike times
    [trial_psth, bins] = hist(ith_trial_times, bins);
    % accumulated the binned times across trials
    binned_trials(i, :) = trial_psth;
end

% compute the psth
psth = sum(binned_trials,1) ./ num_trials ./ params.bin_size; % make units spikes/sec.

% plot the psth
if params.plot_hist
    plot_axes = set_up_fig_or_axes(params.foa);
    plot(bins, psth, 'color', params.hist_color);
    xlabel('seconds')
    ylabel('spike rate');
end
