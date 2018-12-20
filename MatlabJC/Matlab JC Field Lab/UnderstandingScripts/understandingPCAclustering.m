% script to help understand clustering with PCA

% JC 3/7/2016

% make clusters, each row is a cell, each column a different param value

A(1:5,:) = repmat([0,0,0,0,0],[5,1]) ;
A(6:10,:) = repmat([5,5,5,5,5],[5,1]) ;
A(11:15,:) = repmat([50,50,50,50,50],[5,1]) ;
A(16:20,:) = repmat([50,50,50,0,0],[5,1]) ;

A = A+normrnd(0,5,size(A)) ;

% distance matrix
for a=1:size(A,1) ; % for each cell
    for b=1:size(A,1) ;
        D(a,b) = sqrt(sum((A(a,:)-A(b,:)).^2)) ;
    end
end

% PCA on distance metric
[EigVect,EigVal]= eigs(D) ;
    
plot(EigVect(:,1))

