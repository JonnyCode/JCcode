function I = mutualInformation(M, Ps)
%M is the follwingMatrix
%entries are P(r,s)
%rows are different responses
%colums are different stimuli
%
%Ps is a vector of the probability of each stimulus

M(M==0) = eps;
M = M./sum(sum(M));
[Nr, Ns] = size(M);

Pr = sum(M,2);
Pr = Pr./sum(Pr);

I = 0;
for s=1:Ns
    for r=1:Nr
        I = I + M(r,s)*mylog2(M(r,s)/(Pr(r)*Ps(s)));
    end
end

if isnan(I)
    keyboard
end







