function movWhite = movieWhitener(mov,SpatialWhiteFlag,TemporalWhiteFlag)

% JC 10/7/2016

% this function will spatially or temporally whiten a movie, mov(x,y,t)

movWhite = nan(size(mov)) ; % prep matrix
xl = size(mov,1) ;
yl = size(mov,2) ;
lf = size(mov,3) ;
sn = xl * yl ; % number of stixels

% spatial whiten without changing temporal power
if SpatialWhiteFlag ;
    
    % calc mean amplitudes for each orientation and frequency
    movfft = nan(size(size(mov))) ;
    movfft_amp = zeros(size(mov,1),size(mov,2)) ;
    
    for f=1:lf ; % for each frame
        movfft(:,:,f) = fft2(mov(:,:,f)) ; % get frame fft2
        movfft_amp = movfft_amp + abs(movfft(:,:,f))/lf ; % average amplitudes
    end
    
    % average amplitudes across orientations
    movfft_amp_shift = fftshift(movfft_amp) ; % shift so that low freq are at center (dc is center point now)
    movfft_amp_OrientMean = movfft_amp_shift ; % previous matrix as default
    
    x_mat=repmat([1:xl]-ceil(xl/2),yl,1) ;  % x coodinates from center
    y_mat=repmat([1:yl]'-ceil(yl/2),1,xl) ;  % y coodinates
    Dist_mat = sqrt((x_mat-1).^2+(y_mat-1).^2) ; % distance of each frequency from center
    Dist_mat_rnd = round(Dist_mat) ; % estimate bins  
    if xl>=yl ;
        Dist_bins = Dist_mat(1,1:ceil(xl/2)) ;
    else
        Dist_bins = Dist_mat(1:ceil(yl/2),1) ;
    end
    for d=1:length(Dist_bins) ; % for each bin
        movfft_amp_OrientMean(Dist_mat_rnd==Dist_bins(d)) = mean(movfft_amp_shift(Dist_mat_rnd==Dist_bins(d))) ; % average over orientations
    end    
    
    ampMean = mean(movfft_amp_OrientMean(:)) ; % average amplitude across frequencies
    
    % apply a scaling dependent on spatial frequency, but idependent of orientation and frame
    Filter = (1./ifftshift(movfft_amp_OrientMean))*ampMean ; % one filter for entire movie (needed to undo shift first, hence ifftshift)
    for f=1:lf ; % for each frame
        movWhite_fft(:,:,f) = movfft(:,:,f).*Filter  ; % multiple by the Filter 
        movWhite(:,:,f) = real(ifft2(movWhite_fft(:,:,f))) ; 
    end
end
 
% temporal whiten without changing spatial power (is this possible?)
if TemporalWhiteFlag ;
    movfft_amp = zeros(size(mov)) ;
    s=1 ; % stixel count
    for x=1:xl ; % for each x
        for y=1:yl ; % for each y
            movfft(x,y,:) = fft2(mov(x,y,:)) ;
            movfft_amp = movfft_amp+abs(movfft(s,:))/sn ; % average amplitudes
            s=s+1 ;
        end
    end
    ampMean = mean(movfft_amp) ; %   
    Filter = (1./movfft_amp)*ampMean ; % one filter for entire movie
    s=1 ; % stixel count
    for x=1:xl ; % for each x
        for y=1:yl ; % for each y
            movWhite_fft(x,y,:) = movfft(x,y,:).*Filter  ; % multiple by the Filter 
            movWhite(x,y,:) = real(ifft2(movWhite_fft(x,y,:))) ; 
            s=s+1 ;
        end
    end
end
    
    