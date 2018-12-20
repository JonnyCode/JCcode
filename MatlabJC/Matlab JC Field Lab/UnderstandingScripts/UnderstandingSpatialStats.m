% understanding spatial frequency metrics manipulations

% JC 10/5/2016
Spatial_N = 100

% 0 degree sine wave
swf = .01 ;
swA = 1 ;
swp = 0 ;
swo = 0 ;
sw = swA*sind([1:Spatial_N]*swf*360+(swp*360)) + swo ;

sw_mat = repmat(sw,[100,1]) ;
imagesc(sw_mat) % see sine wave
imagesc(abs(fft2(sw_mat))) % see fft amplitudes for a sinewave oriented in one direction
[ps2,ps1,ps1x] = PowerSpectrumFinder2D(sw_mat,1,true); % find power spectrum
imagesc(ps2) ; % see power spectrum
plot(ps1x,ps1) ; % see power spectrum averaged across orientations

% 90 degree sine wave
swf = .2 ;
swA = 1 ;
swp = 0 ;
swo = 0 ;
sw = swA*sind([1:Spatial_N]*swf*360+(swp*360)) + swo ;

sw_mat = repmat(sw',[1,100]) ;
imagesc(sw_mat)

imagesc(abs(fft2(sw_mat))) % see a power spectrum for a sinewave oriented in one direction

% 45 degree sine wave
swf = .2 ;
swA = 1 ;
swp = 0 ;
swo = 1 ;
sw = swA*sind([1:Spatial_N*3]*swf*360+(swp*360)) + swo ;

sw_mat=nan(Spatial_N,Spatial_N) ;
for a=[1:Spatial_N] ;
    for b=0:a-1 ;
        sw_mat(a-b,b+1)=sw(a) ; % distance along the diagonal
    end
end
for a=[1:Spatial_N-1] ;
    for b=0:Spatial_N-a ;
        sw_mat(Spatial_N-b,a+1+b)=sw(a+Spatial_N) ;
    end
end
imagesc(sw_mat) 

imagesc(abs(fft2(sw_mat))) % see a power spectrum for a sinewave oriented in one direction

% look at fft2 of image
temp = load('clown.mat');
im_mat = temp.X ;
im_mat = im_mat-mean(im_mat(:)) ;
imagesc(im_mat)

[ps2,ps1,ps1x] = PowerSpectrumFinder2D(im_mat,1,true); % find power spectrum
imagesc(ps2)  % see power spectrum
plot(ps1x,ps1) ; % see power spectrum averaged across orientations

% gaussian
Gcov = [10,0;0,10] ;
[grid_x,grid_y] = meshgrid([1:Spatial_N],[1:Spatial_N]) ; 
temp = mvnpdf([grid_x(:),grid_y(:)], 50, Gcov) ;
G = reshape(temp, Spatial_N, Spatial_N) ;

[ps2,ps1,ps1x] = PowerSpectrumFinder2D(G,1,true); % find power spectrum
imagesc(ps2)  % see power spectrum
plot(ps1x,ps1) ; % see power spectrum averaged across orientations

% make a fft2 
fft2fake = repmat(.001,Spatial_N)+repmat(.001,Spatial_N)*i ; % create a fake fft2

fft2fake(1,1) = 2+3i ; % choose an orientation in the relative x,y indixes
sw_mat = real(ifft2(fft2fake)) ;
imagesc(sw_mat) ;

