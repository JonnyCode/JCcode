function [Generator,Psth, SNR] = LNmodel_JC(movie_path,datarun, varargin)

movie_path = '/Volumes/lab/acquisition/movie-xml/BW-15-1-0.48-11111-40x40-60.35.xml' ;

% parameters
p = inputParser;

p.addParameter('psthBinTime', 0.1) ; % sec
p.addParameter('nl_bin_num', 20) ; % sec
p.addParameter('displayFrameRate', 60.35) ; % hz
p.addParameter('movieRefresh', 2) ;
p.addParameter('cell_indices','all') ; % datarun indicies of cells to analyze

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
        
        


