function ForIgor = MovieRasterPlots(DataBlock, DB, Params)

% this function will make a raster plot for every repeated movie stimulus
% block that can be mapped over from a BWN stim block

% JC 10/17/2016 

% parameters
Color_list = {'k','r','k','r','k','r','k','r','k','r','k','r'} ; % order of colors for each 
saveFigPath = '/Users/jcafaro/Documents/AnalysisFigures/NatStimRasters/Db2' ;
%stim_time = [10:20] ; % time of stimulus to display - 'all' if you want the whole thing

% binary white data run
% load Master data
dataRunMaster = load_data(DataBlock(DB).BwPath{1}) ;
dataRunMaster = load_neurons(dataRunMaster) ;
dataRunMaster = load_ei(dataRunMaster, 'all') ;
dataRunMaster = load_params(dataRunMaster,'cell_type_depth', 5) ;
dataRunMaster = load_sta(dataRunMaster) ; % only necessary to get trf_time

% load receptive fields and filter
marks_params.thresh = 4.5;
dataRunMaster = get_sta_summaries(dataRunMaster, 'all','marks_params', marks_params);
 
filt_params.radius = 0.75;
dataRunMaster = get_rfs_filtered(dataRunMaster, 'all', 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');

% save spatial and temporal rf so can delete dataRunMaster later
for a=1:length(dataRunMaster.stas.filt_rfs) ;
    srf{a} = dataRunMaster.stas.filt_rfs{a}*dataRunMaster.stas.polarities{a} ; % give the correct pix polarity 
end
trf = dataRunMaster.stas.time_courses ; 

dataRunMaster.stimulus.monitor_refresh = 60.35 ; % the monitor_refresh rate is not accurate
trf_time = -[0:dataRunMaster.stas.depth-1]*dataRunMaster.stimulus.interval/dataRunMaster.stimulus.monitor_refresh ; 

% load and map each movie block data set
for db=1:length(DataBlock(DB).MovieRepPath) ; % for each movie data block
    
    dataRun{db} = load_data(DataBlock(DB).MovieRepPath{db}) ;
    dataRun{db} = load_neurons(dataRun{db}) ;
    dataRun{db} = load_ei(dataRun{db}, 'all') ;
    
    cell_list_map{db} = map_ei(dataRunMaster, dataRun{db}) ;
end

% plot raster
if ~exist('stim_time', 'var') ;
    stim_time =[0, max(DataBlock(DB).MovieRepFrameNum)/60] ;
end

for cells = 1:length(dataRunMaster.spikes) ; % for each cell

    figure(1)
    clf
    for db=1:length(DataBlock(DB).MovieRepPath) ; % for each movie data block
        cell_id = cell_list_map{db}(cells) ;
        if ~isempty(cell_id{1}) ;
            cell_i = get_cell_indices(dataRun{db},cell_id{1}) ;

            trigNum = ceil(DataBlock(DB).MovieRepFrameNum(db)/100) ; % number of triggers per repeat
            trigs = dataRun{db}.triggers(1:trigNum:end) ; % 
            repeatNum(db) = floor(length(dataRun{db}.triggers)/trigNum); % show how many repeats this is

            get_raster(dataRun{db}.spikes{cell_i}, trigs,'tic_color',Color_list{db},...
                'foa',-1, 'first_tic',sum(repeatNum(1:db))-repeatNum(db)+1,'axis_range', [stim_time 0 sum(repeatNum)]) ;
            
            drawnow

        end
    end
    title(['Cell id',num2str(dataRunMaster.cell_ids(cells))])
    print(gcf, '-djpeg', [saveFigPath,'Cell',num2str(cells)])
end

