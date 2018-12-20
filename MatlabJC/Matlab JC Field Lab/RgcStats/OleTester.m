function [DirectionEstimate,DirectionEstimateMag] = OleTester(OleWeights,ResponseMat,varargin)

% this function will test the wieghts from an Optimal linear estimator for a vector quanity.
% JC 8/30/17 (adapted from 'OLE_Angles' script from Joel Zylberberg)

% input
% OleWeights is a matrix (NumCells,2) with an X,Y weight for each cell and possibly one extra 'cell' for constant offset  
% ResponseMat is a Matrix with responses from different neurons in each column and trials in rows.

% output
% DirectionEstimate is a column vector of estimated stimulus direction in degrees 

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

% add fake 'cell' if OleWeights include fake cell weight
if size(ResponseMat,2)+1==size(OleWeights,1) ;
    ResponseMat = [ResponseMat, ones(size(ResponseMat,1),1)] ;
end

% calculate cartisian estimate 
DirEstXY = ResponseMat*OleWeights ; 
        
% polar estimate
DirectionEstimate = atan2d(DirEstXY(:,2), DirEstXY(:,1)) ; % four quadrant inverse tangent
DirectionEstimate(DirectionEstimate<0) = DirectionEstimate(DirectionEstimate<0)+360 ; % no negatives
DirectionEstimateMag = sqrt(DirEstXY(:,2).^2 + DirEstXY(:,1).^2) ;
