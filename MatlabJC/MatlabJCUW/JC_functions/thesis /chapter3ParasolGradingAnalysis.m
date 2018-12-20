function ForIgor = chapter3ParasolGradingAnalysis(Input,A)

% thesis chapter 3 figures of exc and inh in parasol during grading

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);
 
SpatialStimParams = hdf5load(['~/Data/primate/',Input(A).cellname,'_spatial.h5']) ;

epochsExc = str2num(Input(A).Exc) ;
epochsInh = str2num(Input(A).Inh) ;
epochsCA = str2num(Input(A).CA) ;

for a = 1:length(epochsExc) ; % for each g epoch
    [dataExc(a,:),error] = ITCReadEpoch(epochsExc(a), 0, fp) ;    % get data
    
    [SIexc, error] = ITCGetSamplingInterval(epochsExc(a), fp); % get sampling interval
    
    [MonitorExc(a,:), error] = ITCReadEpoch(epochsExc(a), 1 , fp) ;    % get monitor data
    
    StrucString = ['params_epoch',num2str(epochsExc(a))] ;
    Struct = SpatialStimParams.(StrucString) ;
   
    if strcmp('SpatialGradingStimulus', Struct.stimClass) ;

        preFramesExc(a) = Struct.spatial_prepts ;
        stimFramesExc(a) = Struct.spatial_stimpts ;
        postFramesExc(a) = Struct.spatial_postpts ;

        BarWidthExc(a) = Struct.BarWidth ;
        CycleFramesExc(a) = Struct.CycleFrames ;
        PhaseExc(a) = Struct.Phase ;
        AmplitudeExc(a) = Struct.amplitude ;
        MeanLevelExc(a) = Struct.spatial_meanLevel ;
        SpotRadiusExc(a) = Struct.spotRadius ;
        SquareWaveSpaceExc(a) = Struct.squareWaveSpace ;
        SquareWaveTimeExc(a) = Struct.squareWaveTime ;
        sumWavesExc(a) = Struct.sumWaves ;

        SpatialParamsExc(a,:) = [preFramesExc(a),stimFramesExc(a),postFramesExc(a),BarWidthExc(a),CycleFramesExc(a),...
            PhaseExc(a),AmplitudeExc(a),MeanLevelExc(a),SpotRadiusExc(a),SquareWaveSpaceExc(a),SquareWaveTimeExc(a),sumWavesExc(a)] ;
    else
        errmsg('wrong stim')
    end

end

for a = 1:length(epochsInh) ; % for each g epoch
    [dataInh(a,:),error] = ITCReadEpoch(epochsInh(a), 0, fp) ;
 
    [SIinh, error] = ITCGetSamplingInterval(epochsInh(a), fp);
 
    [MonitorInh(a,:), error] = ITCReadEpoch(epochsInh(a), 1 , fp) ;    % get monitor data
    
    StrucString = ['params_epoch',num2str(epochsInh(a))] ;
    Struct = SpatialStimParams.(StrucString) ;
   
    if strcmp('SpatialGradingStimulus', Struct.stimClass) ;

        preFramesInh(a) = Struct.spatial_prepts ;
        stimFramesInh(a) = Struct.spatial_stimpts ;
        postFramesInh(a) = Struct.spatial_postpts ;

        BarWidthInh(a) = Struct.BarWidth ;
        CycleFramesInh(a) = Struct.CycleFrames ;
        PhaseInh(a) = Struct.Phase ;
        AmplitudeInh(a) = Struct.amplitude ;
        MeanLevelInh(a) = Struct.spatial_meanLevel ;
        SpotRadiusInh(a) = Struct.spotRadius ;
        SquareWaveSpaceInh(a) = Struct.squareWaveSpace ;
        SquareWaveTimeInh(a) = Struct.squareWaveTime ;
        sumWavesInh(a) = Struct.sumWaves ;

        SpatialParamsInh(a,:) = [preFramesInh(a),stimFramesInh(a),postFramesInh(a),BarWidthInh(a),CycleFramesInh(a),...
            PhaseInh(a),AmplitudeInh(a),MeanLevelInh(a),SpotRadiusInh(a),SquareWaveSpaceInh(a),SquareWaveTimeInh(a),sumWavesInh(a)] ;
    else
        errmsg('wrong stim')
    end
    
end

SI = SIexc *1e-6;
if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end
time = [1:length(dataExc)]*SI ;

FrameStartPnts = MonitorAnalysis2([MonitorExc; MonitorInh], length(epochsExc)+length(epochsInh), SI) ;
for a = 1:length(FrameStartPnts) ; 
    framerate_mean(a) = 60/mean(diff(FrameStartPnts{a}(2:end))) ; % frames/pnt
end
frameRate = mean(framerate_mean) ;

prePntsExc = floor(preFramesExc/frameRate) ;
postPntsExc = floor(postFramesExc/frameRate) ;
stimPntsExc = floor(stimFramesExc/frameRate) ;

prePntsInh = floor(preFramesInh/frameRate) ;
postPntsInh = floor(postFramesInh/frameRate) ;
stimPntsInh = floor(stimFramesInh/frameRate) ;


% low pass filter to remove electrical crap
samplerate = 1/SI ;
dataExc_lpf = lowPassFilter(dataExc, samplerate, 5000) ; %(signal,samplerate,cutoff frequ (hz))
dataInh_lpf = lowPassFilter(dataInh, samplerate, 5000) ; 

% conductances
for a = 1:length(epochsExc) ;
    gExc = (dataExc - mean(dataExc(a,:)))/-61 ;
end
for a = 1:length(epochsInh) ;
    gInh = (dataInh - mean(dataInh(a,:)))/61 ;
end

gExc = gExc - min(gExc(:)) ;
gInh = gInh - min(gInh(:)) ;

UniqueStim = unique([SpatialParamsExc;SpatialParamsInh],'rows') ;

for a=1:size(UniqueStim,2) ;
    
    UiExc = find(sum(repmat(UniqueStim(a,:),size(SpatialParamsExc,1),1)-SpatialParamsExc,2)==0) ;
    UiInh = find(sum(repmat(UniqueStim(a,:),size(SpatialParamsInh,1),1)-SpatialParamsInh,2)==0) ;
    
    Gexc{a} = gExc(UiExc,:) ;
    Ginh{a} = gInh(UiInh,:) ;

end

for a=1:size(UniqueStim,2) ;
    if ~isempty(Gexc{a}) ;
        Gexc_mean(a,:) = mean(Gexc{a}) ;
    end
    if ~isempty(Ginh{a}) ;
        Ginh_mean(a,:) = mean(Ginh{a}) ;
    end 
    
    if ~isempty(Gexc{a}) & ~isempty(Ginh{a}) ;
        cc(a,:) = xcov(Gexc_mean(a,:),Ginh_mean(a,:),'coeff') ;
    end
end


for a=1:size(UniqueStim,2) ;
    figure
    subplot(4,1,1:3)
    if ~isempty(Gexc{a})
        plot(time,mean(Gexc{a},1),'b')
        hold on
    end
    if ~isempty(Ginh{a})
        plot(time,mean(Ginh{a},1),'r')
    end 
    
    subplot(4,1,4)
    if ~isempty(Gexc{a}) & ~isempty(Ginh{a}) ;
        plot(cc(a,:))
    end
    
end


% For Igor

if A==239 ;
    identifier = ['GexcExample',num2str(A)] ;
    ForIgor.(identifier) = Gexc{12}(1,:);
    
    identifier = ['GinhExample',num2str(A)] ;
    ForIgor.(identifier) = Ginh{12}(1,:);
end
    
    


