function Gen = LinModel_NatStim(moviePath, dataRun, CellId, varargin)

% This function will calculate a Generator signal from the dataRun.sta for
% a natural stimulus movie

% parameters
p = inputParser;

% JCafaro 11/28/16
Color_list = {'r','b','g','c','y','k'} ;

% parse inputs
p = inputParser;
p.addParamValue('GeneratorBinTime', [], @isnumeric) ; % (sec) time bin of generator signal, defaults to movie frame rate (interval*monitor frame rate)
p.addParamValue('MovieFrameInterval', 1, @isnumeric) ; % interval of movie (num of frames displayed for each unique movie frame)
p.addParamValue('MovieFramesPerSecond', 60.35, @isnumeric) ; % frame rate of display
p.addParamValue('MovieStixelWidth', 1, @isnumeric) ; % number of pixels per stixel of the movie
p.addParamValue('MovieX_StartEnd', [0, 800], @isnumeric) ; % number of pixels per stixel of the movie
p.addParamValue('MovieY_StartEnd', [0, 600], @isnumeric) ; % number of pixels per stixel of the movie
p.addParameter('cell_indices','all') ; % datarun indicies of cells to analyze

p.parse(varargin{:});
params = p.Results;

% parameters
if isempty(params.PsthTimeDuration);
    params.PsthTimeDuration = min(diff(PsthStartTimes)) ;
end

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

for c=1:numCells ; % NOT SET FOR c>1 yet!!!!!!!!
    
    % convolve with spatial receptive field % CAUTION RED PIXELS ONLY
    for f=1:movie_frames ; % for each movie frame
        temp = sum(mov(:,:,1,f).*datarun.stas.rfs{params.cell_indices(c)}) ;
        prj_s(f) = sum(temp) ;
    end

    % convolve with temporal receptive field
    tc = [zeros(1,1+length(datarun.stas.time_courses{params.cell_indices(c)})),...
        fliplr(datarun.stas.time_courses{params.cell_indices(c)}')] ;
    prj = conv(prj_s,tc,'same') ;
    
    % upsample generator signal
    movie_si = 1/movieFrameRate ;
    movie_time = [movie_si:movie_si:datarun.duration] ;
    
    si = 1/datarun.sampling_rate ;
    time = [si:si:datarun.duration] ;
    prj_int = interp1(movie_time, prj, time) ;

% get psth acctual data
    spikeTrain = zeros(1,length(time)) ;
    spikePnts = floor(datarun.spikes{params.cell_indices(c)}*datarun.sampling_rate) ;
    spikeTrain(spikePnts) = 1 ;
    
    psth_binSize = params.psthBinTime*datarun.sampling_rate ;
    psth = smooth(spikeTrain,psth_binSize)*psth_binSize/params.psthBinTime ;

% nonlinearity
    prj_int_range = range(prj_int) ;
    prj_bins = [min(prj_int):prj_int_range/params.nl_bin_num:max(prj_int)] ;
    for b = 1:length(prj_bins)-1 ;
        NL(b) = mean(psth(prj_int>=prj_bins(b) & prj_int<prj_bins(b+1))) ;
    end
end
        