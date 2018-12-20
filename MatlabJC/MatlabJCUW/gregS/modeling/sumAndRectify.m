function y = sumAndRectify(x,th,T) 
%T is transformer
y = sum(x);
y = y^2;
if y<th, y = 0; end
    