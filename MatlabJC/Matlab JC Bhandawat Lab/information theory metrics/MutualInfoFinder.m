function [Im,Sentropy] = MutualInfoFinder(S,R,Rstd,Rdelta,Srange,Sprob) 

% JC 4/4/14
% this function will take a set of tuning curves (defined by S,R, and
% Rstd) and calculate the mutual information of R about S (and vice-versa).  
% Rdelta is a constant that sets the precision of R that is trusted (i.e. significant digits of R)  
% Srange is a vector that sets the bins overwhich to inerpolate R 
% Sprob is a vector that sets the probability of values in Srange.  
% S,R,Rstd, should be in arrays for each tuning curve (e.g. S{1},R{1},Rstd{1} should all be from the same tuning curve).

% measure from Dyan and Abbot book page 130
% Im = integrateOverDs(ds*(integrateOverDr((p(s)p(r|s)log2(p(r|s)/p(r)))*dr)) 

% make sure Sprob is pdf and find Sdelta
Sprob = Sprob/sum(Sprob) ;
Sdelta = Srange(2)-Srange(1) ; 

% make R and Rstd percise within bounds
for a=1:length(S) ; % for tuning curve
    R{a} = round(R{a}/Rdelta)*Rdelta ;
    Rstd{a} = round(Rstd{a}/Rdelta)*Rdelta ;
end

% interpolate across R+/-Rstd
for a=1:length(S) ; % for tuning curve
    RstdInterp(a,:) = interp1(S{a},Rstd{a},Srange,'nearest','extrap') ; % nearest neighbor interp Rstd
    Rinterp(a,:) = interp1(S{a},R{a},Srange,'linear','extrap') ; % linear interp R
end
Rinterp(Rinterp<0) = 0 ; % rectify interpolated responses

% calculate p(R|S)
probRbins = [min([Rinterp(:)-RstdInterp(:);0]):Rdelta:max(Rinterp(:)+RstdInterp(:))] ; % the probability bin values of R  
for a=1:length(Srange) ; % for each S bin   
    numObsRgivenS(a,:) = zeros(1,length(probRbins)) ; % start with empty vector 
    for b = 1:length(S) ; % for each tuning curve
        numObsRgivenS(a,:) = numObsRgivenS(a,:) + exp(-((probRbins-Rinterp(b,a)).^2)/(2*RstdInterp(b,a)^2)) ; % gaussian assumption
    end
    probRgivenS(a,:)= numObsRgivenS(a,:)/sum(numObsRgivenS(a,:)) ; %make into pdf  
end

% calculate p(r)=p(s)*p(r|s)
Rprob = probRgivenS'*Sprob' ; % sum the probabilities across values of S and wieght by the probability you will see S at all
Rprob = (Rprob/sum(Rprob))' ; % make pdf

% calculate Im matrix = (p(s)p(r|s)log2(p(r|s)/p(r))) 
for a=1:length(Srange) ; % for each S bin  
    logProb(a,:) = log2(probRgivenS(a,:)./Rprob) ; % the log probability 
    ImMat(a,:) = Sprob(a)*probRgivenS(a,:).*logProb(a,:) ; % the mutual information matrix for each R and S
end

% integral over Im matrix to get mutual information (Im)
ImMat_sumOverR = sum(ImMat*Rdelta,2) ; % integrate over R
Im = sum(ImMat_sumOverR*Sdelta) ; % integrate over S 
    
% calculate stim entropy (Im<=stim entropy)
Sentropy = -sum(Sprob.*log2(Sprob)*Sdelta) ;


    