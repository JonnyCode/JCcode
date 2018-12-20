function [InputBins1,InputBins2,Nonlinearity,Nonlinearity_sem] = NonLinFilterFinder2(Input1,Input2,Output,binsize1,binsize2) ;

% this function will creat an input output curve for two inputs
% JC 8/12/10

% set bins
extrabin1 = rem(range(Input1(:)),binsize1)>0 ;
extrabin2 = rem(range(Input2(:)),binsize2)>0 ;

InputBins1 = min(Input1(:)):binsize1:max(Input1(:))+(binsize1*extrabin1) ;
InputBins2 = min(Input2(:)):binsize2:max(Input2(:))+(binsize2*extrabin2) ;

% make lookup table where data exists
for a = 1:length(InputBins1)-1 ; %for each possible bin
    for b = 1:length(InputBins2)-1 ;
        Outputdata = Output(Input1>=InputBins1(a) & Input1<InputBins1(a+1) & Input2>=InputBins2(b) & Input2<InputBins2(b+1)) ; % find the data within the combo
        Nonlinearity(a,b) = mean(Outputdata) ;
        if length(Outputdata)>1 ;
            Nonlinearity_sem(a,b) = std(Outputdata)/sqrt(length(Outputdata)) ;
        else
            Nonlinearity_sem(a,b) = nan ;
        end
    end
end

% center bins
InputBins1 = InputBins1(1)+binsize1/2:binsize:InputBins1(end) ; 
InputBins2 = InputBins2(1)+binsize2/2:binsize:InputBins2(end) ;

% find nans and interpolate to fill in lookup table where data does not exist
rowInterp = nans(size(Nonlinearity)) ;
for a = 1:length(InputBins1); % for each row 
    nani = find(isnan(Nonlinearity(a,:))) ; % nans 
    if ~isempty(nani) ;
        numi = setdiff([1:length(Nonlinearity(a,:))],nani) ; % non nans 
        rowInterp(a,nani) = interp1(InputBins2(numi),Nonlinearity(a,numi),InputBins2(nani),'linear','extrap') ; % interpolate over nans in each row (this is an interpolation between columns)
    end
end

columnInterp = nans(size(Nonlinearity)) ;
for a = 1:length(InputBins2); % for each column
    nani = find(isnan(Nonlinearity(:,a))) ; % nans 
    if ~isempty(nani) ;
        numi = setdiff([1:length(Nonlinearity(:,a))],nani) ; % non nans 
        columnInterp(nani,a) = interp1(InputBins1(numi),Nonlinearity(numi,a),InputBins1(nani),'linear','extrap') ; % interpolate over nans in each column
    end
end
    
TwoDInterpolation = (rowInterp+columnInterp)/2 ; % avearge over interpolations

Nonlinearity(isnan(Nonlinearity)) = TwoDInterpolation(isnan(Nonlinearity)) ; % fill in nonlinearity