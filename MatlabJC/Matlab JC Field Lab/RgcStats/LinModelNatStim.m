function [GenSig,GenSigTime] = LinModelNatStim(MoviePath, dataRun, CellIndex, varargin)

% This function will calculate a Generator signal from the dataRun.sta for
% a natural stimulus movie

% JCafaro 12/1/16

% parse inputs
p = inputParser;
p.addParamValue('GeneratorBinTime', [], @isnumeric) ; % (sec) time bin of generator signal, defaults to movie frame rate (interval*monitor frame rate)
p.addParamValue('MovieFrameInterval', 1, @isnumeric) ; % interval of movie (num of frames displayed for each unique movie frame)
p.addParamValue('MovieFramesPerSecond', 60.35, @isnumeric) ; % frame rate of display
p.addParamValue('MovieStixelWidth', 1, @isnumeric) ; % number of pixels per stixel of the movie
p.addParamValue('MovieXstart', 0, @isnumeric) ; % start of x position 
p.addParamValue('MovieYstart', 0, @isnumeric) ; % start of y position
p.addParamValue('separateSrfTrfFlag', true, @islogical) ; % false, use full sta; true: construct from timecourse and Srf
p.addParamValue('staPixDimensions', [800,600], @isnumeric) ; % pixel size of BWN used to make sta
p.addParamValue('staFrameInterval', 1, @isnumeric) ; % frame interval of BWN used to make sta
p.addParamValue('Sta',[], @isnumeric) ; % sta (defualt will be calculated from dataRun.stas)

p.parse(varargin{:});
params = p.Results;

% default parameters
MovieTimeStep = params.MovieFrameInterval/params.MovieFramesPerSecond ;
if isempty(params.GeneratorBinTime);
    params.GeneratorBinTime = MovieTimeStep ;
end

% load movie 
Temp = load(MoviePath) ;% load movie
MovieMat=nan(size(Temp.mov,2),size(Temp.mov,1),size(Temp.mov,3)) ; %prep Movie matrix
for t=1:size(Temp.mov,3) ;
    MovieMat(:,:,t) = Temp.mov(:,:,t)' ; % the movie is loaded as transpose of displayed movie
end
clear Temp ;
NumFrames = size(MovieMat,3) ; % number of movie frames

% get sta
if isempty(params.Sta) ; % if no sta was provided
    
    staXnum = size(dataRun.stas.stas{CellIndex},2) ;
    staYnum = size(dataRun.stas.stas{CellIndex},1) ;
    staFnum = size(dataRun.stas.stas{CellIndex},4) ;

    if params.separateSrfTrfFlag ;
        Sta = nan(staYnum,staXnum,staFnum) ;
        SrfNorm = dataRun.stas.rfs{CellIndex}/norm(dataRun.stas.rfs{CellIndex}) ;
        TrfNorm = dataRun.stas.time_courses{CellIndex}/norm(dataRun.stas.time_courses{CellIndex}) ;
        for f=1:staFnum ; % for each frame
            Sta(:,:,f) = SrfNorm * TrfNorm(f) ;
        end
    else
        Sta = dataRun.stas.stas{CellIndex} ;
    end
else
    staXnum = size(Sta,2) ;
    staYnum = size(Sta,1) ;
    staFnum = size(Sta,3) ;
end
    
    
% interpolate sta to match stimulus space/time sample rate
StaStixWidth = floor(params.staPixDimensions(1)/staXnum) ;
StaStixLength = floor(params.staPixDimensions(2)/staYnum) ;

StaX = [.5:staXnum]*StaStixWidth ;
StaY = [.5:staYnum]*StaStixLength ;
StaTime = [0:staFnum-1]*params.staFrameInterval/params.MovieFramesPerSecond ;

FilterXNum = size(MovieMat,2) ;
FilterYNum = size(MovieMat,1) ;

FilterTime = [0:MovieTimeStep:StaTime(end)] ;
FilterX  = [params.MovieXstart+.5:FilterXNum]*params.MovieStixelWidth ;
FilterY  = [params.MovieYstart+.5:FilterYNum]*params.MovieStixelWidth ;

FilterTimeNum = length(FilterTime) ;

% spatial interpolation
staSpaceInterp = nan(FilterYNum, FilterXNum, staFnum) ; % prep
for f=1:staFnum ; % for each frame of sta
    staSpaceInterp(:,:,f)=interp2(StaX,StaY,Sta(:,:,f),FilterX,FilterY','linear',0) ;
end

% temporal interpolation
Filter = nan(FilterYNum,FilterXNum,FilterTimeNum) ;
for x=1:FilterXNum ; % for each X
    for y=1:FilterYNum ; % for each Y
        Filter(y,x,:) = interp1(StaTime,squeeze(staSpaceInterp(y,x,:)),FilterTime) ;
    end
end

% convolve movie with filter
GeneratorSignal = nan(1,NumFrames) ; % prep
for f=1:NumFrames-FilterTimeNum ;
    MovieBlock = MovieMat(:,:,f:f+FilterTimeNum-1) ;
    GeneratorSignal(f+FilterTimeNum) = sum(MovieBlock(:).*Filter(:)) ;
end

GenSigTime = [params.GeneratorBinTime:params.GeneratorBinTime:NumFrames*MovieTimeStep] ;

if params.GeneratorBinTime==MovieTimeStep ;
    GenSig = GeneratorSignal ;
else
    MovieTime = [0:MovieTimeStep:NumFrames*MovieTimeStep] ;
    for f=1:length(GenSigTime)-1 ;
        GenSig(f) = mean(GeneratorSignal(GenSigTime(f)<=MovieTime<GenSigTime(f+1))) ;
    end
    GenSig(length(GenSigTime))=mean(GeneratorSignal(GenSigTime(length(GenSigTime))<=MovieTime)) ;
end
    
 
        