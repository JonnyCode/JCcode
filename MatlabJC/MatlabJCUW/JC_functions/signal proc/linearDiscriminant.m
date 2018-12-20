function pCorrect = linearDiscriminant(covCA,meanM,x,y) 

% 'covCA' is a cell array with a set of covariance matrices. 'meanM' should be a
% 2 column matrix with each row defining the x,y means of each gaussian. 
% 'x' and 'y' are vectors with values at which to evaluate.
%'pCorrect' will indicate the percent correct from an ideal linear discriminant.

% JC 10/18/11

numGauss = length(covCA) ; % number of 2d gaussian distributions to be made
Correct = 0 ;
Total = 0 ;

for a=1:length(x) ;
    for b=1:length(y) ;
        xy = [x(a),y(b)] ; % vector
        
        for c=1:numGauss ; % for each 2d distribtution
           %mVd{c}(a,b) = (1/((2*pi)*det(cov{c})^.5))^(-.5*(xy'-mean(c,:)')'*inv(cov{c})*(xy'-mean(c,:)')); % something is WRONG with this equation
            mVd{c}(a,b) = mvnpdf(xy,meanM(c,:),covCA{c}) ;
            temp(c) = mVd{c}(a,b) ;
        end
            Correct = Correct + max(temp) ; % those parts of the distribution that are correct
            Total = Total + sum(temp) ; % total parts
    end
    
end

pCorrect = Correct/Total ; % percent correct

%figure
% for a=1:numGauss ;
%     subplot(1,numGauss,a)
%     imagesc(mVd{a})
% end

