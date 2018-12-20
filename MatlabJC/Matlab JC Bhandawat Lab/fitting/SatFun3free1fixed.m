function fit = SatFun3free1fixed(coef,Input) ;

% adapted from "SatFun4free"
% fized 2 variables variables

x=Input.x ;
MinY = Input.MinY ;

fit = (coef(3)-MinY)*((x.^coef(1))./(x.^coef(1) +coef(2).^coef(1)))+MinY ;