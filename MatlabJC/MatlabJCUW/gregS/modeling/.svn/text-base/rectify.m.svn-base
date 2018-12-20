function y = rectify(x,th,T) 
%T is transformer
y = T.revert(x);
y(x<th) = 0;
y = T.invert(y);
    