function [pdfMean,CrossTerms_pdfMean,pdfX] = RspProbFinder(DirectionVect,ResponseMat,varargin)

% this function will find pdfs (prob of each direction) for each cell and crossproducts

% input
% DirectionVect is a column vector with direction of stimulus in degrees for each trial.
% ResponseMat is a Matrix with responses from different neurons in each column and trials in rows.
% 'OQE' will use a 'Optimal quadradic estimator' instead of OLE 

% output
% OleWeights provides a set of 'optimal' linear wieghts for each cells that best ('least squares') 
% predicts the DirectionVect, last set is fake cell with stim idependent response.  If 'OQE' than weights 
% will be ordered [W1,...Wn,W12,W13,Wqn,Wfake]
% OlePolarWeights provides wieghts as vectors [direction,magnitude]



if ~isempty(varargin) ; % if there is a varargin
    if strcmp(varargin{1},'OQE') ; % if runing an OQE
        NumCells = size(ResponseMat,2) ; % number cells
        NumTrials = size(ResponseMat,1) ; % number of trials

        J = triu(ones(NumCells),1); % binary array of which elements of product matrix are "not same cell twice"
        ResponseCrossProds = nans(NumTrials,sum(J(:))) ; %prep mat

        for t=1:NumTrials ; % for each trial
            Temp = ResponseMat(t,:)'*ResponseMat(t,:) ; % all products 
            ResponseCrossProds(t,:) = Temp(J==1) ; % only non-same cell cross products
        end

        ResponseMat = [ResponseMat,ResponseCrossProds] ;
    end
end

pdfX = [0:45:180] ;
Dirs = unique(DirectionVect) ;

for cells = 1:NumCells ; % for each cell
    histY(cells,:) = hist(DirectionVect(ResponseMat(:,cells)>0),Dirs) ; % hist
    Temp = PolarVectorAddition([Dirs,histY(cells,:)']); % mean direction
    pdfXCntr = Temp(1) ;
    for d=1:length(Dirs) ; % for each direction
        deltaX(d) = acuteAngle(Dirs(d),pdfXCntr)+rand(1) ; % distance from mean (keeps two distances being exactly the same)
    end
    [dX,di] = sort(deltaX) ;
    pdfTemp = histY(cells,:)/sum(histY(cells,:)) ; %pdf
    pdf(cells,:) = interp1(deltaX(di),pdfTemp(di),pdfX,'linear','extrap') ;
end

pdfMean = mean(pdf,1) ; % mean across cells

for ct = 1:size(ResponseCrossProds,2) ; % for each cross term
    CrossTerms_histY(ct,:) = hist(DirectionVect(ResponseCrossProds(:,ct)>0),Dirs) ; % hist
    Temp = PolarVectorAddition([Dirs,CrossTerms_histY(ct,:)']); % mean direction
    pdfXCntr = Temp(1) ;
    for d=1:length(Dirs) ; % for each direction
        deltaX(d) = acuteAngle(Dirs(d),pdfXCntr)+rand(1) ;  % distance from mean
    end
    pdfTemp = CrossTerms_histY(ct,:)/sum(CrossTerms_histY(ct,:)) ; %pdf
    CrossTerms_pdf(ct,:) = interp1(deltaX,pdfTemp,pdfX,'linear','extrap') ; 
end
    
CrossTerms_pdfMean = mean(CrossTerms_pdf,1) ; % mean across cells
    





