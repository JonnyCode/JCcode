function MovieConverter(moviePath,NewMovieName,varargin) 

% gets movies into mat file and raw movie formats and will also adjust mean 
% and std to specifications 

% JC 11/28/16 

% movie_path must end in .avi or .mp4 (eg '/Volumes/lab/Documents/Movies/squirrel_video.mp4')
% mat file and raw movie will be saved in same folder as original movie

p = inputParser;
p.addParamValue('RawMovieMean', [], @isnumeric) ; % mean to make raw movie ([1-255])
p.addParamValue('RawMovieStd', [], @isnumeric) ; % std to make raw movie 

p.parse(varargin{:});
params = p.Results;

% getting movies into mat file
Slashi = strfind(moviePath,'/') ;
base_folder = moviePath(1:Slashi(end)) ;
matPath = [base_folder,NewMovieName,'.mat'] ;
RawMoviePath = [base_folder,NewMovieName,'.RawMovie'] ;

mvObj = VideoReader(moviePath) ;
mv = mvObj.read ;
%movie_player(mv) ;

% getting movies into raw.movie format

mv = squeeze(mv(:,:,1,:)) ; % reduce rgb to Intensity

% change to double
mv = double(mv) ;

numFrames = size(mv,3) % display the number of frames

mov = nan(size(mv,2),size(mv,1),size(mv,3)) ; % prep mat for transpose
for f=1:numFrames ; % for each frame
    mov(:,:,f) = mv(:,:,f)' ; % transpose
end

if ~isempty(params.RawMovieMean) ; % adjust mean
    mov = mov-mean(mov(:))+params.RawMovieMean ; 
end

if ~isempty(params.RawMovieStd) ; % adjust std
    mov = (mov/std(mov(:)))*params.RawMovieStd ; 
end

mov(mov>255) = 255 ;
mov(mov<0) = 0 ;

save(matPath,'mov') ;
    
% write movie into raw movie format
stixel = 1 ;

write_movie(matPath,RawMoviePath,stixel)



% option 2
%     mvObj = VideoReader(movie_path) ;
%     mvWidth = mvObj.Width ;
%     mvHeight = mvObj.Height ;
% 
%     mv = struct('cdata',zeros(mvHeight,mvWidth,3,'uint8'),'colormap',[]);
% 
%     k = 1;
%     %while hasFrame(xyloObj)
%     while k<10 ;
%         mv(k).cdata = readFrame(mvObj);
%         k = k+1;
%     end
% 
%     hf = figure;
%     movie(hf,mv,1,mvObj.FrameRate);

    
% option 2
%     mov = nan(size(mv,1),size(mv,2),size(mv,4)) ;
%     for f=1:size(mv,4) ; 
%         mmov(:,:,f) = mean(mv(:,:,:,f),3) ;
%     end

