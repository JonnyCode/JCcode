% to change matlab movie into avi

myVideo = VideoWriter('myfile2.avi'); % save in cd
open(myVideo);
writeVideo(myVideo, M); % Movie file is M
close(myVideo)