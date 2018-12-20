function ForIgor = NatStimMappedBw(DataBlock, DB, Params) 

% will test and LN model on natural stimulus movie
% JC

FramesPerTrig = 100 ;

% load data
dataRun = load_data(DataBlock(DB).MovieRepPath{Params.MovieNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun, 'all') ;

% load data master
dataRunMaster = load_data(DataBlock(DB).BwPath{1}) ;
dataRunMaster = load_neurons(dataRunMaster) ;
dataRunMaster = load_ei(dataRunMaster, 'all') ;
dataRunMaster = load_params(dataRunMaster,'cell_type_depth', 5) ;
dataRunMaster = load_sta(dataRunMaster) ;

for c=1:length(dataRunMaster.cell_types) ; % for each cell type
    Temp{c} = dataRunMaster.cell_types{c}.name ;
end
CellTypeNames = unique(Temp) ;

marks_params.thresh = 4.5;
dataRunMaster = get_sta_summaries(dataRunMaster, 'all','marks_params', marks_params);

map = map_ei(dataRunMaster,dataRun) ;

MasteriMapped = cell(1,length(CellTypeNames)) ;
CellIdsArray = cell(1,length(CellTypeNames)) ;
CelliArray = cell(1,length(CellTypeNames)) ;

for CellType = 1:length(CellTypeNames) ; % for each cell type
    Masteri{CellType} = get_cell_indices(dataRunMaster,CellTypeNames{CellType}) ; % cell index in master
    for c=1:length(Masteri{CellType}) ; % for each cell
        if ~isempty(map{Masteri{CellType}(c)}) ;
            MasteriMapped{CellType} = [MasteriMapped{CellType},Masteri{CellType}(c)] ;
            CellIdsArray{CellType} = [CellIdsArray{CellType},map(Masteri{CellType}(c))] ;
            temp = map(Masteri{CellType}(c)) ;
            CelliArray{CellType} = [CelliArray{CellType},get_cell_indices(dataRun,temp{1})] ;
        end
    end
end

MoviePath = DataBlock(DB).MoviePath{Params.MovieNum} ;

% calculate movie triggers
trigNum = ceil(DataBlock(DB).MovieRepFrameNum(Params.MovieNum)/FramesPerTrig) ; % number of triggers per repeat
PsthStartTimes = dataRun.triggers(1:trigNum:end) ; % begining of each movie 

stixWidth = DataBlock(DB).MovieStixWidth(Params.MovieNum) ; % stixel width
MovieInterval = DataBlock(DB).MovieFrameInterval(Params.MovieNum) ; % frame interval

% for ei oriented movie of firing rates
load(DataBlock(DB).TformEiPath) ;

StimPsthMovie = StimPsthMovieMaker(dataRun, CellIdsArray, Tform, ...
    MoviePath, PsthStartTimes, 'MovieFrameInterval', MovieInterval,...
    'MovieStixelWidth', stixWidth) ;

% test LN model

TestCellType = 27 ;
TestCell = 2;

[GenSig,GenSigTime] = LinModelNatStim(MoviePath, dataRunMaster, MasteriMapped{TestCellType}(TestCell),...
    'MovieFrameInterval', MovieInterval,'MovieStixelWidth', stixWidth, 'separateSrfTrfFlag', true ) ;

[psth,psthTime]=get_smooth_psth(dataRun.spikes{CelliArray{TestCellType}(TestCell)},PsthStartTimes,'stop',min(diff(PsthStartTimes)),'bin_size',0.0331) ;

GenSigTimeEndPnt = find(GenSigTime<=psthTime(end),1,'last') ; % last time point in psth
psthInterp = interp1(psthTime,psth,GenSigTime(1:GenSigTimeEndPnt),'linear','extrap') ;

LnBinNumber = 20 ;
[InputBins,Nonlinearity,Nonlinearity_sem] = NonLinFilterFinder(GenSig(1:GenSigTimeEndPnt),psthInterp,range(GenSig)/LnBinNumber) ;

LnPrediction = interp1(InputBins,Nonlinearity,GenSig,'linear','extrap') ;

% exp fit for NL
% StrPnt=20 ; % to avoid nans in exp fit
% expParams = fit(GenSig(StrPnt:GenSigTimeEndPnt)',psthInterp(StrPnt:end)','exp1') ; % fit exp
% LnPrediction = expParams.a*exp(expParams.b*GenSig(1:GenSigTimeEndPnt)) ;


figure
plot(psthTime,psth,'k')
hold on
plot(GenSigTime(1:GenSigTimeEndPnt),LnPrediction(1:GenSigTimeEndPnt),'r')
plot(GenSigTime(1:GenSigTimeEndPnt),GenSig(1:GenSigTimeEndPnt)*max(psth)/max(GenSig(1:GenSigTimeEndPnt)),'g')

figure
plot(InputBins,Nonlinearity)

figure
norm_rf = norm_image(dataRunMaster.stas.rfs{MasteriMapped{TestCellType}(TestCell)});
imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 

figure
plot(dataRunMaster.stas.time_courses{MasteriMapped{TestCellType}(TestCell)})
