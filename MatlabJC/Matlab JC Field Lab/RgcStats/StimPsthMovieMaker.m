function StimPsthMovie = StimPsthMovieMaker(dataRun, cellIdsArray, Ei2MonTform, ...
    MoviePath, PsthStartTimes, varargin)

% this function will plot the location of the cells center of mass and a
% circle proportional to the cells firing rate

% JCafaro 11/28/16
Color_list = {'r','b','g','c','y','k'} ;

% parse inputs
p = inputParser;
p.addParamValue('MonitorBounds', [0, 800, 0, 600], @isnumeric) ; % pixels (region of interst in monitor space)
p.addParamValue('PsthBinTime', 0.1, @isnumeric) ; % (sec) time of psth bin
p.addParamValue('PsthNormalizeFlag', false, @islogical) ; % normalize psth by max of firing rate
p.addParamValue('SynchronyLineFlag', false, @islogical); % plot lines between cells when they fire at the same time
p.addParamValue('PsthTimeDuration', [],@isnumeric) ; % time duration of psth (defaults to shortest time between start times)
p.addParamValue('StimPsthMovieTimeStep', 0.016, @isnumeric) ; % time step of the compiled psth and stim movie
p.addParamValue('MovieFrameInterval', 1, @isnumeric) ; % interval of movie (num of frames displayed for each unique movie frame)
p.addParamValue('MovieFramesPerSecond', 60.35, @isnumeric) ; % frame rate of display
p.addParamValue('MovieStixelWidth', 1, @isnumeric) ; % number of pixels per stixel of the movie
p.addParamValue('MovieXstart', 0, @isnumeric) ; % start of x position 
p.addParamValue('MovieYstart', 0, @isnumeric) ; % start of y position
p.addParamValue('numElectrodeLayers', 2, @isnumeric) ; % number of electrode layers around single electrode used to get ei center of mass 
p.addParamValue('SpikeDisplayScale', .1, @isnumeric) ; % ~ proportional to number of pixels displayed per spike 

p.parse(varargin{:});
params = p.Results;
    
% parameters
if isempty(params.PsthTimeDuration);
    params.PsthTimeDuration = min(diff(PsthStartTimes)) ;
end

numTypes = length(cellIdsArray) ; % for each cell type

for Type = 1:numTypes ;
    numCells(Type) = length(cellIdsArray{Type}) ;
end

% get psth
for Type = 1:numTypes ; % for each type of cell
    for Cells = 1:numCells(Type) ; % for each cell of that type
        ci = get_cell_indices(dataRun,cellIdsArray{Type}(Cells)) ; % cell index
        [psth{Type}(Cells,:), PsthTime] = get_psth(dataRun.spikes{ci}, PsthStartTimes ,...
                'stop',params.PsthTimeDuration,'bin_size',params.PsthBinTime) ;

        if params.PsthNormalizeFlag ;
            psth{Type}(Cells,:)= psth{Type}(Cells,:)/max(psth{Type}(Cells,:)) ;
        end    
    end
end
   
% get ei center locations
for Type = 1:numTypes ; % for each type
    for Cells = 1:numCells(Type) ; % for each cell of that type
        EiCom{Type}(Cells,:) = get_ei_com(dataRun, cellIdsArray{Type}(Cells), params.numElectrodeLayers) ;
    end
end
    
% load movie 
Temp = load(MoviePath) ;% load movie
MovieMat=nan(size(Temp.mov,2),size(Temp.mov,1),size(Temp.mov,3)) ; %prep Movie matrix
for t=1:size(Temp.mov,3) ;
    MovieMat(:,:,t) = Temp.mov(:,:,t)' ; % the movie is loaded as transpose of displayed movie
end
clear Temp ;
NumFrames = size(MovieMat,3) ; % number of movie frames

% Time 
CompMovieTime = [0:params.StimPsthMovieTimeStep:params.PsthTimeDuration] ; % time points of movie/psth composite
MovieTime = [0:NumFrames-1]*(params.MovieFrameInterval/params.MovieFramesPerSecond) ; % time points of movie

% Movie scaling and offset
MovieX = [params.MovieXstart+1,(params.MovieXstart+size(MovieMat,2))*params.MovieStixelWidth] ;
MovieY = [params.MovieYstart+1,(params.MovieYstart+size(MovieMat,1))*params.MovieStixelWidth] ;

figure
ColorRandFlag = true ; % plot different color for each cell regardless of "Type"
RandColori = randi(length(Color_list),1,numCells) ;
for t=1:length(CompMovieTime) ; % for each time point in composite image
    PsthPnt = interp1(PsthTime,[1:length(PsthTime)],CompMovieTime(t),'nearest','extrap') ;
    MoviePnt = interp1(MovieTime,[1:NumFrames], CompMovieTime(t),'nearest','extrap') ;
    
    colormap gray
    imagesc(MovieX,MovieY,MovieMat(:,:,MoviePnt),[0,255]) ; % plot movie frame in grey scale
    hold on
    
    for Type = 1:numTypes ; % for each type
        for Cells = 1:numCells(Type) ; % for each cell of that type

            Resp = params.SpikeDisplayScale * psth{Type}(Cells,PsthPnt) ; % scale size of response
            
            ctrMon = tformfwd(Ei2MonTform,EiCom{Type}(Cells,:)) ; % ei--> monitor coordinates

            if ~ColorRandFlag ;
                plot(ctrMon(1),ctrMon(2),'+','Color',Color_list{Type}) 
            else
                plot(ctrMon(1),ctrMon(2),'+','Color',Color_list{RandColori(Cells)}) 
            end
            hold on
            if Resp>0 ;
                if ~ColorRandFlag ;
                    plot(ctrMon(1),ctrMon(2),'o','Color',Color_list{Type},'MarkerSize',Resp) 
                else
                    plot(ctrMon(1),ctrMon(2),'o','Color',Color_list{RandColori(Cells)},'MarkerSize',Resp) 
                end
            end
            
            if params.SynchronyLineFlag ; %
                for Cells2 = 1:numCells(Type) ; % for each cell of that type
                    if Cells2~=Cells ;
                        Resp = params.SpikeDisplayScale/10 *sqrt(psth{Type}(Cells,PsthPnt)*psth{Type}(Cells2,PsthPnt)) ; % geometric mean of responses
                        if Resp>0 ;
                            ctrMon2 = tformfwd(Ei2MonTform,EiCom{Type}(Cells2,:)) ; % ei--> monitor coordinates
                            plot([ctrMon(1),ctrMon2(1)],[ctrMon(2),ctrMon2(2)],'-','Color',Color_list{Type},'LineWidth',Resp)
                        end
                    end
                end
            end   
        end
    end
    axis(params.MonitorBounds)
    StimPsthMovie(t)=getframe ;
    hold off
end


        
