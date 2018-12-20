function ForIgor = PNadaptationAnalysisEpop(ForIgor) 

% this function will analyze data across multiple cell outputs from
% PNadaptationAnalysisE
% JC 10/28/12

% average across delta spike data and fit a sigmoid to data

field = fieldnames(ForIgor) ;

b=0 ; % contrations used in each cell
for a=1:length(field) ;
    if strncmp('LogConcentrationcell',field{a},20) ;
        b = b+1 ;
        concentrations{b} = ForIgor.(field{a}) ;
        
    end
end
ConcentrationsUnique = unique(cell2mat(concentrations)) ;

deltaSRcontrol = nan(length(concentrations),length(ConcentrationsUnique)) ;
deltaSRexp = nan(length(concentrations),length(ConcentrationsUnique)) ;
deltaSRwash = nan(length(concentrations),length(ConcentrationsUnique)) ;
deltaSRcontrolStd = nan(length(concentrations),length(ConcentrationsUnique)) ;
deltaSRexpStd = nan(length(concentrations),length(ConcentrationsUnique)) ;
deltaSRwashStd = nan(length(concentrations),length(ConcentrationsUnique)) ;


b=0 ; % arrange data into consistent concentration dependent matricies
for a=1:length(field) ;
    if strncmp('DeltaSRmean1cell',field{a},16) ;
        b = b+1 ;
        [c,ri] = intersect(ConcentrationsUnique',concentrations{b}') ;
        deltaSRcontrol(b,ri) = ForIgor.(field{a}) ;
    end
end

b=0 ;
for a=1:length(field) ;
    if strncmp('DeltaSRmean2cell',field{a},16) ;
        b = b+1 ;
        [c,ri] = intersect(ConcentrationsUnique',concentrations{b}') ;
        deltaSRexp(b,ri) = ForIgor.(field{a}) ;
    end
end

b=0 ;
for a=1:length(field) ;
    if strncmp('DeltaSRmean3cell',field{a},16) ;
        b = b+1 ;
        [c,ri] = intersect(ConcentrationsUnique',concentrations{b}') ;
        deltaSRwash(b,ri) = ForIgor.(field{a}) ;
    end
end

b=0 ;
for a=1:length(field) ;
    if strncmp('DeltaSRstd1cell',field{a},15) ;
        b = b+1 ;
        [c,ri] = intersect(ConcentrationsUnique',concentrations{b}') ;
        deltaSRcontrolStd(b,ri) = ForIgor.(field{a}) ;
    end
end

b=0 ;
for a=1:length(field) ;
    if strncmp('DeltaSRstd2cell',field{a},15) ;
        b = b+1 ;
        [c,ri] = intersect(ConcentrationsUnique',concentrations{b}') ;
        deltaSRexpStd(b,ri) = ForIgor.(field{a}) ;
    end
end

b=0 ;
for a=1:length(field) ;
    if strncmp('DeltaSRstd3cell',field{a},15) ;
        b = b+1 ;
        [c,ri] = intersect(ConcentrationsUnique',concentrations{b}') ;
        deltaSRwashStd(b,ri) = ForIgor.(field{a}) ;
    end
end

% weighted average for all concentrations tested
deltaSRcontrolMean = nanmean(deltaSRcontrol.*(1./deltaSRcontrolStd.^2),1)./nanmean((1./deltaSRcontrolStd.^2),1) ;
deltaSRexpMean = nanmean(deltaSRexp.*(1./deltaSRexpStd.^2),1)./nanmean((1./deltaSRexpStd.^2),1) ;
deltaSRwashMean = nanmean(deltaSRwash.*(1./deltaSRwashStd.^2),1)./nanmean((1./deltaSRwashStd.^2),1) ;

% fit with sigmoid (fiting with sigmoid is prob not good idea until
% saturation points at both low and high end are found

    
    




