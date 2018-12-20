function [RFspat,RFtemp,RFtemp_X] = get_RF_JC(movie_path,spikes_times, varargin)

% getting 

movie_path = '/Volumes/lab/acquisition/movie-xml/BW-15-1-0.48-11111-40x40-60.35.xml' ;

% parameters
p = inputParser;

p.addParameter('displayFrameRate', 60.35) ; % hz
p.addParameter('movieRefresh', 2) ;
p.addParameter('Movie_lag_time',0) ; % (sec) time movie lags behind spike times ((+)--> movie lags, (-)--> spikes lag)  

p.parse(varargin{:});

params = p.Results;

% from parameters
movieFrameRate = params.displayFrameRate / params.movieRefresh;
movie_frames = floor(movieFrameRate * datarun.duration);

if isnumeric(params.cell_indices) ;
    numCells = length(params.cell_indices) ;
else
    numCells = length(datarun.spikes) ;
end

% get movie
mov = get_movie(movie_path, datarun.triggers, movie_frames);
