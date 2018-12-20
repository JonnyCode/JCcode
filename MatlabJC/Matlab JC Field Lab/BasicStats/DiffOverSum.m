function deltaIndex = DiffOverSum(A,B) 

% this function will calculate A-B/(abs(A)+abs(B))

deltaIndex = (A-B)./(abs(A)+abs(B)) ;