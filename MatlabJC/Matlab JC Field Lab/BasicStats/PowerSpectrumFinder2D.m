function [PowerSpectrum2D,VectorSpectrum,VectorSpectrumX] = ...
    PowerSpectrumFinder2D(Signal, SampleRate, VectorSpectrumFlag)

% JC 10/7/2016 
% This function will take a set of frames, Signal(x,y,t), and a SampleRate in (Stix/um),
% and calculate a 2d spatial power spectum and a estimate collapsed across orientations
% called 'VectorSpectrum' if VectorSpectrumFlag=true.

xl = size(Signal,1) ;
yl = size(Signal,2) ;
fn = size(Signal,3) ; % frame number
sn = xl*yl ; % number of stixel per frame

% 2d power spectrum (high frequencies are in the center, low frequencies in the corners, DC is upper left corner)
for f=1:fn ;
    Sfft2 = fft2(Signal(:,:,f)) ; % fourier transform
    PsFrame(:,:,f) = (abs(Sfft2)/sn).^2 ; % power spectrum normalized by number of stixel (not by variance)
end
PowerSpectrum2D = mean(PsFrame,3) ; % average across all frames

% Colapse across orientations (Avoiding points in fft2 beyond nyquist in longest orientation)
if VectorSpectrumFlag ;
    Niquist = SampleRate/2 ;
    points = max(xl,yl) ;
    VectorSpectrumX = [0:2*Niquist/points:Niquist] ;
    lVsX = length(VectorSpectrumX) ;
    
    x_mat=repmat([1:xl],yl,1);  % x coodinates
    y_mat=repmat([1:yl]',1,xl);  % y coodinates
    Dist_mat = sqrt((x_mat-1).^2+(y_mat-1).^2) ; % distance of each stixel from center
    Dist_mat_rnd = round(Dist_mat) ; % estimate bins
    
    if xl>=yl ;
        Dist_bins = Dist_mat(1,1:lVsX) ;
    else
        Dist_bins = Dist_mat(1:lVsX,1) ;
    end
    for d=1:length(Dist_bins) ; % for each bin
        VectorSpectrum(d) = mean(PowerSpectrum2D(Dist_mat_rnd==Dist_bins(d))) ; % average over orientations
    end    
else
   VectorSpectrum = nan ;
   VectorSpectrumX = nan ;
end

