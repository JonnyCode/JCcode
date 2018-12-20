function fit = HockeyStickFree2(coef,x)

% coef 1 = inflection point
% coef 2 = slope 

fit=nan(1,length(x)) ;
for a=1:length(x) ;
    if x(a)<coef(1) ;
        fit(a) = 0 ;
    else
        fit(a:end) = coef(2)*(x(a:end)-coef(1)) ;
    end
end

