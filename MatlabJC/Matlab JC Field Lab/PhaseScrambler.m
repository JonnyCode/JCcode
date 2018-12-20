function movScrambled = movPhaseScrambler(mov)

% JC 10/6/2016 
% this function will take a movie(x,y,t) of pixel intesity and phase
% scramble so that the spatial and temporal power spectum of the movie are
% unchanged

xl = size(mov,1) ; % length x
yl = size(mov,2) ; % length y
fl = size(mov,3) ; % length frames

% load random image to get consisten phase noise that has appropriate symetry
rng(1) ; % set random seed so it is consistent
RandomImage = rand(xl,yl) ;

RandomImage_fft2 = fft2(RandomImage) ; % fft
RandomImage_phase = angle(RandomImage_fft2) ; % random phase noise that has correct symetry for fft2

% load movie and shift each frequency the same on each frame
for f = 1:fl ; % for each frame
    Fft2 = fft2(mov(:,:,f)) ;
    Fft2_amp = abs(Fft2) ; % magnitude
    Fft2_phase = angle(Fft2) ; % phases
    Fft2_phase_delta = Fft2_phase + RandomImage_phase ; % add random image phase shift
    Fft2_scrambled = Fft2_amp.*exp(1i*Fft2_phase_delta) ; 
    movScrambled(:,:,f) = real(ifft2(Fft2_scrambled)) ;
end


    