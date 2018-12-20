function [OleWeights, OlePolarWeights] = OleFinder(DirectionVect,ResponseMat,varargin)

% this function will find wieghts for a Optimal linear estimator for a vector quanity.
% JC 8/30/17 (adapted from 'OLE_Angles' script from Joel Zylberberg)

% input
% DirectionVect is a column vector with direction of stimulus in degrees for each trial.
% ResponseMat is a Matrix with responses from different neurons in each column and trials in rows.
% 'OQE' will use a 'Optimal quadradic estimator' instead of OLE
% 'Regularize' will penalize big wieghts by 'lambda' eg- OleFinder('OQE','Regularize',10)

% output
% OleWeights provides a set of 'optimal' linear wieghts for each cells that best ('least squares') 
% predicts the DirectionVect, last set is fake cell with stim idependent response.  If 'OQE' than weights 
% will be ordered [W1,...Wn,W12,W13,Wqn,Wfake]
% OlePolarWeights provides wieghts as vectors [direction,magnitude]


% put directionVect into X,Y
DirXY(:,1) = cosd(DirectionVect) ;
DirXY(:,2) = sind(DirectionVect) ;

if ~isempty(varargin) ; % if there is a varargin
    if strcmp(varargin{1},'OQE') ; % if runing an OQE
        NumCells = size(ResponseMat,2) ; % number cells
        NumTrials = size(ResponseMat,1) ; % number of trials

        J = triu(ones(NumCells),1); % binary array of which elements of product matrix are "not same cell twice"
        ResponseCrossProds = nans(NumTrials,sum(J(:))) ; %prep mat

        for t=1:NumTrials ; % for each trial
            Temp = ResponseMat(t,:)'*ResponseMat(t,:) ; % all products 
            %Temp = (repmat(ResponseMat(t,:)',1,NumCells)-repmat(ResponseMat(t,:),NumCells,1))./...
            %(repmat(ResponseMat(t,:)',1,NumCells)+repmat(ResponseMat(t,:),NumCells,1)+10^-8) ; % all ratio quotients (alternative to products)
            ResponseCrossProds(t,:) = Temp(J==1) ; % only non-same cell cross products
        end

        ResponseMat = [ResponseMat,ResponseCrossProds] ;
    end
end

% add 'cell' with constant response so that Ole could fit a constant offset
Resp = [ResponseMat,ones(size(ResponseMat,1),1)] ; 

% find Ole Wieghts
if length(varargin)>1 
    if strcmp(varargin{2},'Regularize') ; % if regularizing (penalize bigger weights) 
        lambda = varargin{3} ; % penalty parameter for large wights (try 0.1-100, 0 should be same as not regularizing)
        OleWeights(:,1) = LassoIteratedRidge(Resp,DirXY(:,1),lambda) ;
        OleWeights(:,2) = LassoIteratedRidge(Resp,DirXY(:,2),lambda) ;
    end
else
    OleWeights = Resp\DirXY ;
end

% Put Ole wieghts into vectors with [degrees, magnitude]
PolarDir = atan2d(OleWeights(:,2), OleWeights(:,1)) ; % four quadrant inverse tangent
PolarDir(PolarDir<0) = PolarDir(PolarDir<0)+360 ; % no negatives
PolarMag = sqrt(sum(OleWeights.^2,2)) ; % magnitude
OlePolarWeights = [PolarDir,PolarMag] ;

% Below - to figure out which OQE cross terms belong to which cells
% CrossTermCelli = [] ;
% for cells = 2:NumCells ;
%     CrossTermCelli = [CrossTermCelli,[repmat(cells,[1,cells-1]);[1:cells-1]]] ; % cellA;cellB in order of crossterms
% end


