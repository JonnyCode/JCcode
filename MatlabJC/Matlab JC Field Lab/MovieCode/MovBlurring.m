function movBlurred = MovBlurring(mov,StdSpatialBlur,StdTemporalBlur)

% JC 10/7/2016
% this function will apply a guassian, defined by Std,
% to blur a movie, mov(x,y,t).  If only leave zero if you don't want to
% blur space or time; blurs space first

movBlurred = nan(size(mov)) ; 

xl = size(mov,1) ;
yl = size(mov,2) ;
fl = size(mov,3) ;

if StdSpatialBlur>0 ; % blur space

    CovMat = [StdSpatialBlur^2,0;0,StdSpatialBlur^2] ;
    [grid_x,grid_y] = meshgrid([-xl:xl],[-yl:yl]) ; 
    temp = mvnpdf([grid_x(:),grid_y(:)], [0,0], CovMat) ;
    Filter = reshape(temp, yl*2+1, xl*2+1)' ;
    
    for f=1:fl ;
        movBlurred(:,:,f) = conv2(mov(:,:,f),Filter,'same') ;
    end
    
else
    movBlurred = mov ;
end
   
if StdTemporalBlur>0 ; % blur space
 
    FilterTd = pdf('norm',[-fl:fl], 0, StdTemporalBlur) ; % filter in the time domain
    
    for x=1:xl ; 
        for y=1:yl ;
            movBlurred(x,y,:) = conv(squeeze(movBlurred(x,y,:)),FilterTd,'same') ; % convolve with filter (not sure why shift is necessary)
        end
    end
    
end
