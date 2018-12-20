
%% examples on how to run spectra

raw_data_path = '/Volumes/backup008/2016-02-18-0/data004';
analysis_path = '/Analysis/gdfield/2016-02-18-0/data004';
movie_path = '/Volumes/lab/acquisition/movie-xml/BW-15-1-0.48-11111-40x40-60.35.xml';

mVision(raw_data_path, analysis_path, '',movie_path, 'all','none','')

mVision(raw_data_path, analysis_path, '',movie_path, [0,0,0,0,0,0,1,1], [0,0,0,0,0,0,1,1],'')

%%
/Applications/MATLAB_R2015a.app/bin/matlab -nodesktop % for terminal after ssh

cd /Volumes/Haydn/Users/circuit/Development/spectra ;
raw_data_path = '/Volumes/data/2016-10-26-0/data009';
analysis_path = '/Volumes/Haydn/Analysis/jcafaro/2016-10-26-0/data009';
movie_path = '/Volumes/lab/acquisition/movie-xml/BW-15-2-0.48-11111-40x40-60.35.xml';

mVision(raw_data_path, analysis_path, '',movie_path, 'all','none','')




