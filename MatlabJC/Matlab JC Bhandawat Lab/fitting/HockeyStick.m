function fit = HockeyStick(coef,x)

% coef 1 = inflection point
% coef 2 = offset
% coef 3 = slope 

fit=nan(1,length(x)) ;
for a=1:length(x) ;
    if x(a)<coef(1) ;
        fit(a) = coef(2) ;
    else
        fit(a:end) = coef(3)*(x(a:end)-coef(1))+coef(2) ;
    end
end

