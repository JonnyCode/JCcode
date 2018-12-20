function [jpsth] = JpsthFinder(X,Y,varargin) 

% will find the joint psth for a set of trials in rows X and Y matricies.
% Adding 'Coef' to last option will normalize the JSPTH to corr coef

% JC 10/24/12

if size(X,1)~=size(Y,1) ;
    error('X and Y must have same number of trials')
end
    
J = zero(size(X,2),size(Y,2)) ;
Xmean = mean(X,1) ; 
Ymean = mean(Y,1) ;
Xstd = std(X,0,1) ;
Ystd = std(Y,0,1) ;

for a=size(X,1) ; % for each trial
    for b=1:size(X,2) ; % for each time point in X
        for c=1:size(Y,2) ; % for each time point in Y
            J(b,c) = Jcov(b,c)+(X(a,b)-Xmean(b))*(Y(a,c)-Ymean(c)) ;
        end
    end
end
J = J/size(X,1) ;

if strcmp(varargin{1},'Coef')
    for b=1:size(X,2) ; % for each time point in X
        for c=1:size(Y,2) ; % for each time point in Y
            J(b,c) = J(b,c)/(Xstd(b)*Ystd(c)) ;
        end
    end
end


