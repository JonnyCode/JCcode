clear all

load -ASCII c10p1.mat

%scatter(c10p1(:,1),c10p1(:,2))

mean_c10p1 = mean(c10p1)
c10p1_adjust(:,1) = c10p1(:,1)-mean_c10p1(1,1)
c10p1_adjust(:,2) = c10p1(:,2)-mean_c10p1(1,2)
scatter(c10p1_adjust(:,1),c10p1_adjust(:,2))
hold on

tau_w = 1
delta=.01
w_trace = []
u = c10p1_adjust
w = rand(1,2) ; 
alpha = 1

for trial = 1:6
w = rand(1,2) ; 
for round=1:100
for i=1:100

v=u(i,:)*w' ;
    
%delta_w = (1/tau_w)*(v*u(i,:)-(alpha*v^2)*w) ;
delta_w = (v*u(i,:))-((v^2)*w) ;

w = w + delta*delta_w ;

w_trace=[w_trace w'] ;

end


end

plot(w_trace(1,:),w_trace(2,:),'r-')
hold on, plot(w_trace(1,:),w_trace(2,:),'o')
end
xlabel('w1')
ylabel('w2')

% figure, plot(diff(w_trace(1,:))) ;
% figure, plot(diff(w_trace(2,:))) ;

cor_matrix = cov(u) ;
[vector,value]=eig(cor_matrix) ;

