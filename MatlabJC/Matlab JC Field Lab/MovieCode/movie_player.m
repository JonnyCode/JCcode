function movie_player(movie) 

movie_frames  = size(movie,3) ;

figure
for f=1:movie_frames ; % for each movie frame
    imagesc(movie(:,:,f))
    pause
end