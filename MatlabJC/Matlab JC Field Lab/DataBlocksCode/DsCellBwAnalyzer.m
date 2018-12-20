function ForIgor = DsCellBwAnalyzer(DataBlock, DB, Params)

numDs = length(DataBlock(DB).BwPath) ;

for ds = 1:numDs % for each data set
    dataRun = load_data(DataBlock(DB).BwPath{ds}) ;
    dataRun = load_neurons(dataRun) ;
    dataRun.stimulus.monitor_refresh = 60.35 ; % the monitor_refresh rate is not accurate
    
    marks_params.thresh = 3 ;
    dataRun = load_sta(dataRun) ;
    dataRun = load_params(dataRun) ;
    dataRun = get_sta_summaries(dataRun, 'all','marks_params', marks_params) ;
    
    
    % get movie
    dataRun = load_java_movie(dataRun) ;
    mv = get_movie(DataBlock(DB).BwMoviePath{ds}, dataRun.triggers, size(dataRun.stimulus.java_movie)) ;

    % get spike triggered variance
    stv = spikeTrigVarBw(dataRun,[1:2],mv) ;
    stv = spikeTrigVarBlur(dataRun,1,mv, 3, 3) ;
    
    trf_time(ds,:) = -[0:dataRun.stas.depth-1]*dataRun.stimulus.interval/dataRun.stimulus.monitor_refresh ; % time axis of temporal receptive field
    PixPerStix = dataRun.stimulus.stixel_height ;
end