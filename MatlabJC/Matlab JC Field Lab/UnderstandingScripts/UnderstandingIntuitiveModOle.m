% imagine two cells that often spike together

q = [1,.8;.8,1] ; % covariance mat
W = [0,0;3,3] ; % weights
R{1} = [3,3]' ; % both cells active
R{2} = [3,0]' ; % only the untuned cell active
R{3} = [0,3]' ; % only the tuned cell is active
alpha = .2 ; 
beta = .2 ;

for a=1:length(R) ;
    ole = R{a}'*W ;
    mod_alpha = -alpha*(q*R{a}.^2)'*W ;
    mod_beta1 = beta*diag((q*R{a})'*(R{a}.*repmat(W(:,1),1,size(R{a},2)))) ;
    mod_beta2 = beta*diag((q*R{a})'*(R{a}.*repmat(W(:,2),1,size(R{a},2)))) ;
    est{a} = ole + mod_alpha + [mod_beta1,mod_beta2] ;
end

