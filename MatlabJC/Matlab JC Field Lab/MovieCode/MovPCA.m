function [movPCs,movPrj,PcVal] = MovPCA(mov,FrameNumber, varargin)

% JC 4/27/2017 
% this function will run PCA on a sets of sequential movie frames
% mov(x,y,t)
% FrameNumber = number of frames in a bin (should be ~ stim autocorr)

p = inputParser;
p.addParamValue('SlidingWindowFlag', true, @islogical); % might avoid the sliding window if ensemble gets too big
p.parse(varargin{:});
params = p.Results;

SlidingWindowFlag = params.SlidingWindowFlag ;

% get ensemble of movie parts 
xl = size(mov,1) ; % x length
yl = size(mov,2) ; % y length
fl = size(mov,3) ; % frame length

if SlidingWindowFlag ; % if sliding window
    mel = fl - FrameNumber + 1 ; % movie ensemble length 
else % if independent 
    mel = floor(fl/FrameNumber) ; % movie ensemble length 
end

MovEnsemble_size = [mel,xl*yl*FrameNumber] ; % size of ensemble when mov is linearized
MovEnsemble = nans(MovEnsemble_size) ; % prep mat

for f=1:mel ; % for each set in the ensemble
    if SlidingWindowFlag ; % if sliding window
        strpnt = f ;    
    else
        strpnt = (f-1)*FrameNumber+1 ; % start point
    end
    Temp = mov(:,:,strpnt:strpnt+FrameNumber-1) ; % movie
    MovEnsemble(f,:) = Temp(:) ; % linearize
end

% run PCA and project
[movPCs,movPrj,PcVal] = pca(MovEnsemble,'Economy',true); % economy deals with the large scale of the cov matrix
PcVal = PcVal/sum(PcVal) ;

% movPCs - PCs are in descending order arranged in columns
% to reconstruct into movie space --> PCn = reshape(movPCs(:,n),[xl,yl,FrameNumber])

% movPrj - each column is a PC, each row is a projection of the stim onto that PC


% figure
% plot(ValFract,'*')
    


