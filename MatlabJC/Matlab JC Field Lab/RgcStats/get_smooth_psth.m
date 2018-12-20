function [psth,psthTime] = get_smooth_psth(spike_times, trigger_times, varargin)

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
% varargin 'tailTime': time at beginging and end of psth that should be
% guessed (by default this is half the bin width).
%

% SET UP OPTIONAL ARGUMENTS
p = inputParser;
p.addParamValue('start', 0, @isnumeric);
p.addParamValue('stop', [], @isnumeric);
p.addParamValue('sample_rate', 20000, @isnumeric);
p.addParamValue('first_trial', 1, @isnumeric)
p.addParamValue('plot_hist', false, @islogical);
p.addParamValue('bin_size', 0.1, @isnumeric);
p.addParamValue('hist_color', [0 0 0]);
p.addParamValue('foa', 0);
p.addParamValue('labels', true, @islogical);
p.addParamValue('tailTime', [], @isnumeric); 

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
psthTime = [params.start:1/params.sample_rate:params.stop] ;
psthTime = round(psthTime*params.sample_rate)/params.sample_rate ; % round within sample rate (not sure why necessary)
SpikeTrain = zeros(num_trials,length(psthTime)) ;

for trg = 1:num_trials ; % for each trigger
    spike_times_tc = spike_times - trigger_times(trg + params.first_trial -1) ; % trigger centered
    spike_times_tc = round(spike_times_tc * params.sample_rate)/params.sample_rate ; % round within sample rate (not sure why necessary) 
    [v,i] = intersect(psthTime, spike_times_tc) ;
    SpikeTrain(trg,i) = 1 ;
end
SpikeTrain_mean = mean(SpikeTrain,1) ; % average spike train
SpikeTrain_mean_smooth = smooth(SpikeTrain_mean,params.bin_size*params.sample_rate) ;

% compute the psth
if isempty(params.tailTime)
    tailPnts = ceil(params.bin_size*params.sample_rate/2) ; % points before full smooth bin is found
else
    tailPnts = ceil(params.tailTime*params.sample_rate) ;
end

psth = SpikeTrain_mean_smooth * params.sample_rate ; % make units spikes/sec.
psth(1:tailPnts) = psth(tailPnts+1) ; % guess tail points
psth(end-tailPnts+1:end) = psth(end-tailPnts-1) ;

% plot the psth
if params.plot_hist
    plot_axes = set_up_fig_or_axes(params.foa);
    plot(bins, psth, 'color', params.hist_color);
    xlabel('seconds')
    ylabel('spike rate');
end
