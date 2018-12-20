
load('/Volumes/lab/MouseCam')

figure
for mv=1:2 ;
    clf
    set(gcf,'name',['movie:',num2str(mv)])
    for a=1:size(MouseCam{mv}.movie,3); 
        imagesc(MouseCam{mv}.movie(:,:,a)); 
        colormap(gray); 
        pause 
    end
end

% for avi movie
obj = VideoReader('/Volumes/lab/Movie5_1stpart.avi');
vid = read(obj);
obj.NumberOfFrames