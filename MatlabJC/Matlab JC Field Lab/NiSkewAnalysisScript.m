% analysis of image skew in natural images


SrfWidth = [20:20:150].^2 ; % (pix) variance of srf
fs = 400 ;

for a=1:length(SrfWidth) ; % for each filter size
    [grid_x,grid_y] = meshgrid([1:fs],[1:fs]) ; 
    temp = mvnpdf([grid_x(:),grid_y(:)], [fs,fs]/2, [SrfWidth(a),0;0,SrfWidth(a)]) ;
    Srf{a} = reshape(temp, fs, fs)/sum(temp(:)) ;
    
    for im = 1:length(NatImages) ;
        lp = conv2(NatImages{im},Srf{a},'same') ; % linear prediction
        lpMedianMinusMean(im,a) = median(lp(:))-mean(lp(:)) ; % median - mean 
    end
end
 
for im = 1:length(NatImages) ;
    imMedianMinusMean(im) = median(NatImages{im}(:))-mean(NatImages{im}(:)) ; % % median - mean 
end

figure
for im = 1:length(NatImages) ; 
    plot([0,sqrt(SrfWidth)],[imMedianMinusMean(im),lpMedianMinusMean(im,:)])
    hold on
end
xlabel('Filter std (pix)')
ylabel('median - mean of lp')
