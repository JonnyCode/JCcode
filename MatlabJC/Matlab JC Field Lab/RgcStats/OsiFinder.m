function [Osi_vl,Osi_ah,Po_ah] = OsiFinder(directions,response) 

% JC 10/3/2016
% This function will find and Orientation Selective Index (Osi) based on vector length
% from Mazurek et al 2014 "Robust quantification..." and based on an ad hoc
% guess from JC.
% the perfered orientation (Po) from the product of the symetric response.

% 'directions' should be in radians and symetric

% sort in ascending order and make row vectors
[d,di]=sort(directions) ;
r=response(di) ;
d = reshape(d,1,length(d)) ;
r = reshape(r,1,length(r)) ;

Osi_vl = abs(sum(r.*exp(2*1i*d))/sum(r)) ; % vector length from Mazurek et al 2014

% ad hoc calculations
PairPnts = length(d)/2 ;
r_pd = r(1:PairPnts).*r(PairPnts+1:end) ;

Osi_ah = max(r_pd)/sum(r_pd) ; % the largest symetric product/sum of symetric responses 

[m,mi]= max(r_pd) ;
Po_ah = d(mi) ;
    
    
    


