function [LinFilters,Constant] = MultiCellLinFilterFinder(StimVector,ResponseMat,FilterLength)

% JC 5/23/18

% This function will find the linear filters to perform optimal linear decoding.  
% This is different than finding the Linear Filter of each cell individually.
% Taken from Warland et al 1997 equation 8.

% Stim vector = stimulus that your training on s(time) in column vector
% ResponseMat is a Matrix with responses from different neurons in each column and stim points in rows.
% FilterLength = number of points in the filter

% make sure Stime vector is column
StimVector = StimVector(:) ;

% params
sl = length(StimVector) ;
fl = FilterLength ;
nc = size(ResponseMat,2) ; % number of cells
rl = sl-fl ; % size of R

% construct R (rows are time points, columns are different cells or same cell with time shifted response)
R = nan(rl,fl) ; % prep mat
R(:,1) = ones(rl,1) ; % first column is constant

cnt = 2 ;
for c=1:nc ; % for each cell
    for f = 1:fl ; % for each point in filter
        R(:,cnt) = ResponseMat(f:f-1+rl,c) ; % shift response by one point
        cnt=cnt+1 ; % add to count
    end
end

% find LinFilters
f = (R'*R)\(R'*StimVector(1:rl)) ; % Warland equation 8

% divide into matrix with each column a cell and each row a filter point
Constant = f(1) ;
LinFilters = reshape(f(2:end),fl,nc) ;
        


