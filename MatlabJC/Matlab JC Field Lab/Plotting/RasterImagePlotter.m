function RasterImagePlotter(SpikeTimes,StartTimes,varargin) 

% plotting rasters from array of spike times 'SpikeTimes' at time points,
% 'StartTimes'.

% JC 2018-11-20

% SET UP OPTIONAL ARGUMENTS
p = inputParser;
p.addParamValue('start', 0, @isnumeric); % time added to StartTime (negative is before start time)
p.addParamValue('stop', [], @isnumeric); % time added to StartTime (should only be positive)
p.addParamValue('TimeStep',0.001, @isnumeric); % (sec) time scale resolution
p.addParamValue('tic_color', [0 0 0]);
p.addParamValue('line_space', 0.2, @isnumeric)
p.addParamValue('tic_thickness', 1);
p.addParamValue('first_trial', 1, @isnumeric)
p.addParamValue('foa', 0)
p.addParamValue('plot', true);
p.addParamValue('labels', true);
p.addParamValue('axis_range', [0 20 0 50]);
p.addParamValue('first_tic', 1, @isnumeric) ;
p.addParamValue('LineTrials',[], @isnumeric) ; % trial numbers for separating lines

% parse inputs
p.parse(varargin{:});
params = p.Results;


%%% function begins here %%%

% figure out the duration of each trial
if isempty(params.stop)
    params.stop = min(diff(StartTimes));
end

% params
Duration = params.stop-params.start ; % (sec)
NumBins = floor(Duration/params.TimeStep) ; % number of bins
NumTrials = length(StartTimes) ;

% round spike and start times to time step resolution
SpikeTimes = round(SpikeTimes/params.TimeStep)*params.TimeStep ;
StartTimes = round(StartTimes/params.TimeStep)*params.TimeStep ;

% make sparse matrix
rasterPlot = sparse(NumTrials,NumBins) ; % sparse mat

for t=1:NumTrials ; % for each trial
    spks = SpikeTimes - StartTimes(t) ; % times relative to start time t
    spks = spks(spks>=params.start & spks<params.stop) ; % just those within the window
    spkPnts = 1+floor((spks-params.start)/params.TimeStep) ; % spike points 
    rasterPlot(t,spkPnts) = 1 ;
end

% plot raster
% figure
imagesc([params.start:params.TimeStep:params.stop],[1:NumTrials],~rasterPlot)
colormap gray
set(gca,'YDir','normal') ; % trial 1 at bottom

hold on % add lines
for t=1:length(params.LineTrials) ; 
    plot([0,NumBins],[params.LineTrials(t),params.LineTrials(t)])
end


