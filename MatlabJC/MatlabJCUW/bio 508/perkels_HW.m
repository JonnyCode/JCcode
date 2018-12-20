% Perkel HW 11/10/07
%Q1
F = 5
K = 2
n = 4
Ca = [.001:.001:100];

A = F*((Ca./(Ca+K)).^n);


plot(Ca,A,'b')
h = gca ;
set(h,'xscale','log')
xlabel ('[Ca] mM')
ylabel ('EPSP amplitude (mV)')
gtext([num2str(F),num2str(K),num2str(n)])
hold on

% Q2

F = 5
K = 2
n = 4
Ca = [.001:.001:1];

A = F*((Ca./(Ca+K)).^n);


plot(Ca,A,'g.-')
h = gca ;
set(h,'xscale','log')
set(h,'yscale','log')
xlabel ('[Ca] mM')
ylabel ('EPSP amplitude (mV)')
gtext([num2str(F),num2str(K),num2str(n)])
hold on

slope = (log10(A(1000))-log10(A(1)))/(log10(Ca(1000))-log10(Ca(1)))

%Q2
Ca = [.4,1,1.8,6]
A = [17.9, 23.5, 26.6, 43.3]
A = [0.1, 7, 16.1, 25]
A = [.9, 4, 5.9, 11.4]
A = [.8, 3.3, 6.6, 11.8]

n = (log10(A(2)) - log10(A(1)))./(log10(Ca(2))-log10(Ca(1)))
F =A(4)

K = ((Ca.*nthroot(F,n))./nthroot(A,n))-Ca


plot(Ca,A) ;
h = gca ;
set(h,'xscale','log')
xlabel ('[Ca] mM')
ylabel ('EPSP amplitude (mV)')


% Q4 problem set 2
n = 5
a = 8
p = .1
count =0

for p = 0:.01:1
    count =count+1
predmean(count) = n*p*a
predvar(count) = n*p*(1-p)*a^2
predCVinv(count) = (predmean(count)^2)/predvar(count)
end

for amp = 0:a:n*a
    count = count+1
prob(count) = (factorial(n)*(p^(amp/a))*((1-p)^(n-(amp/a))))/(factorial(amp/a)*factorial(n-(amp/a))) ;
end

figure,
plot([0:a:n*a],prob)

% q8 ec
%static
n = 100
p = .1
a = 8

% p = normrnd(,sigma_exc,1,length) ;


for trial = 1:100

a = rand(1,n)*16 ;
meana(trial) =mean(a) ;

p = rand(1,n) ;
meanp(trial) = mean(p) ;

transRel = rand(1,n) ;
fire = (p-transRel)>0 ;

amplitude(trial) = sum(fire.*a) ;

end

meanamp = mean(amplitude)
predmean = n*mean(meanp)*mean(meana)

varamp = var(amplitude)
predvar = n*mean(meanp)*(1-mean(meanp))*(mean(meana))^2

