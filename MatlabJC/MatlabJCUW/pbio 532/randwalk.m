%function CV=randwalk(m,n)
ps =[];
mean = [];
var = [];
CV =[];
for i=1:10

p=i*0.1;

%initialize parameters
m=30; %#exc states
n=10; %# inh states
p= 0.1 %probability of going on step up = probabiliy that EPSP arrives
q=1-p; %probability of going on step down = probabiliy that IPSP arrives

%transition matrix
P = diag(q*ones(m+n,1),1) + diag(p*ones(m+n,1),-1); P(1,:) = 0; P(1,1) = 1; P(end,end)=q;

%calculate mean, variance and CV of number of steps to reach threshold with canonical form of P
Q=P(2:m+n+1,2:m+n+1)
I=eye(size(Q));
N=inv(I-Q);

%meanvector=N * ones(m+n,1);
meanvector=sum(N,2)
means(i)=meanvector(m)
%means(i)=sum(N(m,:))

varvector= (2*N-I)* meanvector - meanvector.^2;
vars(i)=varvector(m);

CV(i)= ((vars(i))^0.5)/means(i);
ps(i) = p;
end

plot(ps,CV,'r+:')
figure; plot(CV,'b+:')