function fit = SatFun2free2fixed(coef,Input) ;

% adapted from "SatFun4free"
% fized 2 variables variables

x=Input.x ;
MinY = Input.MinY ;
MaxY = Input.MaxY ;

fit = (MaxY-MinY)*((x.^coef(1))./(x.^coef(1) +coef(2).^coef(1)))+MinY ;