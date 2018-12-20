f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk00399.imc', 'rb', 'ieee-be');
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
colormap(gray);
imagesc(buf');

AxisStd = 100 ;
SmoothNumFrames = 4 ;
Output = ImageJitter(buf, AxisStd,'SmoothNumFrames',SmoothNumFrames) ; 
save('/Users/jcafaro/Desktop/Mov_Image00399Shift100_smooth4.mat','Output','-V7.3') ;

AxisStd = 100 ;
SmoothNumFrames = 1 ;
Output = ImageJitter(buf, AxisStd,'SmoothNumFrames',SmoothNumFrames) ; 
save('/Users/jcafaro/Desktop/Mov_Image00399Shift100_smooth1.mat','Output','-V7.3') ;

AxisStd = 100 ;
SmoothNumFrames = 4 ;
Output = ImageJitter(buf, AxisStd,'SmoothNumFrames',SmoothNumFrames, 'NumFrames', 60*60*10,'AxisShiftMax', 300) ; 
save('/Users/jcafaro/Desktop/Mov_Image00399Shift100_smooth4.mat','Output','-V7.3') ;

% saved for Image_Jitter
f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk00001.iml', 'rb', 'ieee-be'); 
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
VhImage = double(buf) ;
VhImage = VhImage-min(VhImage(:)) ;
VhImage = VhImage*(255/max(VhImage(:))) ; % scale fto fit 0-255
VhImage = VhImage' ;
imagesc(VhImage)
save('/Users/jcafaro/Documents/VhNaturalImages/VhImage_imk00001.mat','VhImage') ; % trees and building

f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk00002.iml', 'rb', 'ieee-be'); 
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
VhImage = double(buf) ;
VhImage = VhImage-min(VhImage(:)) ;
VhImage = VhImage*(255/max(VhImage(:))) ; % scale fto fit 0-255
VhImage = VhImage' ;
colormap gray
imagesc(VhImage,[0,255])
save('/Users/jcafaro/Documents/VhNaturalImages/VhImage_imk00002.mat','VhImage') ; % trees

f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk00408.iml', 'rb', 'ieee-be'); 
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
VhImage = double(buf) ;
VhImage = VhImage-min(VhImage(:)) ;
VhImage = VhImage*(255/max(VhImage(:))) ; % scale fto fit 0-255
VhImage = VhImage' ;
imagesc(VhImage)
save('/Users/jcafaro/Documents/VhNaturalImages/VhImage_imk00408.mat','VhImage') ; % grass and trees

f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk00544.iml', 'rb', 'ieee-be'); 
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
VhImage = double(buf) ;
VhImage = VhImage-min(VhImage(:)) ;
VhImage = VhImage*(255/max(VhImage(:))) ; % scale fto fit 0-255
VhImage = VhImage' ;
imagesc(VhImage)
save('/Users/jcafaro/Documents/VhNaturalImages/VhImage_imk00544.mat','VhImage') ; % water pic

f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk01010.iml', 'rb', 'ieee-be'); 
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
VhImage = double(buf) ;
VhImage = VhImage-min(VhImage(:)) ;
VhImage = VhImage*(255/max(VhImage(:))) ; % scale fto fit 0-255
VhImage = VhImage' ;
imagesc(VhImage)
save('/Users/jcafaro/Documents/VhNaturalImages/VhImage_imk01010.mat','VhImage') ; % too dark

f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk01370.iml', 'rb', 'ieee-be'); 
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
VhImage = double(buf) ;
VhImage = VhImage-min(VhImage(:)) ;
VhImage = VhImage*(255/max(VhImage(:))) ; % scale fto fit 0-255
VhImage = VhImage' ;
imagesc(VhImage)
save('/Users/jcafaro/Documents/VhNaturalImages/VhImage_imk01370.mat','VhImage') ; % tall grass

f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk01600.iml', 'rb', 'ieee-be'); 
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
VhImage = double(buf) ;
VhImage = VhImage-min(VhImage(:)) ;
VhImage = VhImage*(255/max(VhImage(:))) ; % scale fto fit 0-255
VhImage = VhImage' ;
imagesc(VhImage)
save('/Users/jcafaro/Documents/VhNaturalImages/VhImage_imk01600.mat','VhImage') ; % open court and building


f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk01750.iml', 'rb', 'ieee-be'); 
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
VhImage = double(buf) ;
VhImage = VhImage-min(VhImage(:)) ;
VhImage = VhImage*(255/max(VhImage(:))) ; % scale fto fit 0-255
VhImage = VhImage' ;
imagesc(VhImage)
save('/Users/jcafaro/Documents/VhNaturalImages/VhImage_imk01750.mat','VhImage') ; % dirt ground

f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk01923.iml', 'rb', 'ieee-be'); 
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
VhImage = double(buf) ;
VhImage = VhImage-min(VhImage(:)) ;
VhImage = VhImage*(255/max(VhImage(:))) ; % scale fto fit 0-255
VhImage = VhImage' ;
imagesc(VhImage)
save('/Users/jcafaro/Documents/VhNaturalImages/VhImage_imk01923.mat','VhImage') ; % dark tree cover

f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk02016.iml', 'rb', 'ieee-be'); 
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
VhImage = double(buf) ;
VhImage = VhImage-min(VhImage(:)) ;
VhImage = VhImage*(255/max(VhImage(:))) ; % scale fto fit 0-255
VhImage = VhImage' ;
imagesc(VhImage)
save('/Users/jcafaro/Documents/VhNaturalImages/VhImage_imk02016.mat','VhImage') ; % dark tree cover and light tree cover

f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk03002.iml', 'rb', 'ieee-be'); 
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
VhImage = double(buf) ;
VhImage = VhImage-min(VhImage(:)) ;
VhImage = VhImage*(255/max(VhImage(:))) ; % scale fto fit 0-255
VhImage = VhImage' ;
imagesc(VhImage)
save('/Users/jcafaro/Documents/VhNaturalImages/VhImage_imk03002.mat','VhImage') ; % trees

f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk04210.iml', 'rb', 'ieee-be'); 
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
VhImage = double(buf) ;
VhImage = VhImage-min(VhImage(:)) ;
VhImage = VhImage*(255/max(VhImage(:))) ; % scale fto fit 0-255
VhImage = VhImage' ;
imagesc(VhImage)
save('/Users/jcafaro/Documents/VhNaturalImages/VhImage_imk04210.mat','VhImage') ; % too dark to tell

f1 = fopen('/Volumes/Macintosh HD/Users/jcafaro/Downloads/imk04208.iml', 'rb', 'ieee-be'); 
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
VhImage = double(buf) ;
VhImage = VhImage-min(VhImage(:)) ;
VhImage = VhImage*(255/max(VhImage(:))) ; % scale fto fit 0-255
VhImage = VhImage' ;
imagesc(VhImage)
save('/Users/jcafaro/Documents/VhNaturalImages/VhImage_imk04208.mat','VhImage') ; % sand hole (mouse hole?)


