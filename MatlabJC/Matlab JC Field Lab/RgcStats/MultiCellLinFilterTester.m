function StimEstimate = MultiCellLinFilterTester(LinFilters,Constant,ResponseMat)

% this function will test a set of linear filters from
% "MultiCellLinFilterFinder'.  ResponseMat has a different column for each
% cell (in order of LinearFilters).

% JC 5/24/18

% params
sl = size(ResponseMat,1) ; % length of stim est
fl = size(LinFilters,1) ;
nc = size(ResponseMat,2) ; % number of cells
rl = sl-fl ; % size of R

% make R
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

% construct f
f = [Constant;LinFilters(:)] ;

% estimate Stim
StimEstimate = R*f ;
