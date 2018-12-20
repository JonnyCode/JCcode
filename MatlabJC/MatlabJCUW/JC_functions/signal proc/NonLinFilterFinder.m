function [InputBins,Nonlinearity,Nonlinearity_sem] = NonLinFilterFinder(Input,Output,binsize) ;

% this function will create a simple binned input/output relationship curve
% "Input" and "Output" should be the same size arranged in matrix with each
% element in Input corresponding to an element in Output.
% JC 8/12/10

% error messages
if size(Input,1)~=size(Output,1) | size(Input,2)~=size(Output,2) ;
    error('size of input must match size of output')
end

% set bins
extrabin = rem(range(Input(:)),binsize)>0 ;
InputBins = min(Input(:)):binsize:max(Input(:))+(binsize*extrabin) ;

% make lookup table where data exists
for a = 1:length(InputBins)-1 ; %for each possible bin
    Outputdata = Output(Input>=InputBins(a) & Input<InputBins(a+1)) ; % find the data within the combo
    Nonlinearity(a) = mean(Outputdata) ;
    if length(Outputdata)>1 ; 
        Nonlinearity_sem(a) = std(Outputdata)/sqrt(length(Outputdata)) ;  % find sem
    else 
        Nonlinearity_sem(a) = nan ;
    end
end

% center bins
InputBins = InputBins(1)+binsize/2:binsize:InputBins(end) ; 

% find nans and interpolate to fill in lookup table where data does not exist
nani = find(isnan(Nonlinearity)) ; % nans 
if ~isempty(nani) ;
    numi = setdiff([1:length(Nonlinearity)],nani) ; % non nans 
    Nonlinearity(nani) = interp1(InputBins(numi),Nonlinearity(numi),InputBins(nani),'linear','extrap') ; % interpolate over nans in each row (this is an interpolation between columns)
end


    