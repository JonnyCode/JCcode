function Lp = LinDirFilters(Movie, Params) 

% JC 9/29/17

% this function will create a set of linear directional filters, defined by
% Params.  Convolve them with a movie (x,y,t).

% params
Cntrs = Params.centers ; % (pix) X,Y positions of filter centers (Nx2), 
Rad = Params.radi ; % (pix) radius of filter
DirPref = Params.direction_preferences ; % (deg or Nan) direction preference of each filter (Nx1)
Spd = Params.speed ; % (pix/frame) speed preference of filters
SrfFlag = Params.SrfFlag ; % true if it should have a gaussian shape
% Params.Trf - mat of temporal receptive field (time wieghting profile), each row is a filter  

fs = ceil(Rad*2) ; % number of pixels of filter
fl = ceil(fs/Spd) ; % number of filter frames
nF = length(DirPref) ; % number of filters
CycleWidth = fs*2 ;
SrfWidth = fs*4 ; 

DirPref_Unique = unique(DirPref(~isnan(DirPref))) ; % unique directions no nans
if sum(isnan(DirPref))>0 ; % if there is a nan
    DirPref_Unique = [DirPref_Unique,nan] ; % include it
end
    
nft = length(DirPref_Unique) ; % number of filter templates

% spatial receptive field
if SrfFlag ;
    [grid_x,grid_y] = meshgrid([1:fs],[1:fs]) ; 
    temp = mvnpdf([grid_x(:),grid_y(:)], [fs,fs]/2, [SrfWidth,0;0,SrfWidth]) ;
    Srf = reshape(temp, fs, fs)/max(temp) ;
end

% make one filter template for each direction
FltT = cell(1,nft) ;
for d = 1:length(DirPref_Unique) ; % for each unique direction preference
    FltT{d} = nan(fs,fs,fl) ; % prep filter template
    for f=1:fl ; % for each frame
        for x = 1:fs ; % for each x pix coordinate
            for y = 1:fs ; % for each y pix coordinate
                if ~isnan(DirPref_Unique(d)) ; % if its directional
                    FltT{d}(x,y,f) =  cosd((x*cosd(DirPref_Unique(d)) + y*sind(DirPref_Unique(d)) + ((f-1)*Spd))*360/CycleWidth+fs/2)  ;
                    % cos(2dimensional angled line * cyclewidth + phase offset)
                elseif isnan(DirPref_Unique(d)) ; % if non ds
                    FltT{d}(x,y,f) = 1 ; % a square filter
                end 
            end
        end
        if SrfFlag ;
            FltT{d}(:,:,f) = FltT{d}(:,:,f).*Srf ; 
        end
        if isfield(Params,'Trf') ; % temporal receptive field
            FltT{d}(:,:,f) = FltT{d}(:,:,f)*Params.Trf(d,f) ;
        end
    end
end

% make filters with spatial locations
%parpool(4) ; % open parallel pool 

Flt = cell(1,length(DirPref)) ; % prep cell array 
for F = 1:nF ; % for each filter
    Flt{F} = zeros(size(Movie,1),size(Movie,2),fl) ; % prep filter
    if ~isnan(DirPref(F))    
        dpi = find(DirPref_Unique==DirPref(F)) ; % dir pref index
    elseif isnan(DirPref(F)) ;
        dpi = isnan(DirPref_Unique) ;
    end
    
    xpos = [Cntrs(F,1)-floor(fs/2)+1:Cntrs(F,1)+floor(fs/2)] ; % x
    xpos(xpos<1) = 1 ; % in case its off the edge
    xpos(xpos>size(Movie,1)) = size(Movie,1) ;
    
    ypos = [Cntrs(F,2)-floor(fs/2)+1:Cntrs(F,2)+floor(fs/2)] ; % y
    ypos(ypos<1) = 1 ; % in case its off the edge
    ypos(ypos>size(Movie,2)) = size(Movie,2) ;
    
    for f=1:fl ; % for each frame
        Flt{F}(xpos,ypos,f) = FltT{dpi}((xpos>0 & xpos<size(Movie,1)+1),(ypos>0 & ypos<size(Movie,2)+1),f) ;
    end
end
    
for F = 1:nF ; % for each filter
    Flt{F} = Flt{F} - mean(Flt{F}(:)) ;
end

% convolve filter with Movie 
Lp = zeros(nF,size(Movie,3)) ;
for F = 1:length(DirPref) ; % for each filter
    disp(F) ; % for sanity
    for frm = 1:size(Movie,3)-fl ;
        temp = Movie(:,:,frm:frm+fl-1).*Flt{F} ;
        Lp(F,frm+fl-1) = sum(temp(:)) ;
    end
end
       
% figures

% figure % filters
% FltNum = 3 ;
% for f=1:fl ;
%     imagesc(Flt{FltNum}(:,:,f))
%     colorbar
%     pause
% end
%           
% figure % linear predictions
% plot([1:length(Lp)],Lp.^2)
            
            