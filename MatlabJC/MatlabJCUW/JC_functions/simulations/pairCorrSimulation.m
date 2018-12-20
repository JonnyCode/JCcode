% simulating noise correatlions in pairs
for a=1:100 ;
    

lt=10000 ;

% independant waverforms
Ec = normrnd(0,rand,1,lt) ; % common only between cells
Ic = normrnd(0,rand,1,lt) ;

I1 = normrnd(0,rand,1,lt) ; % totally indepedant
E1 = normrnd(0,rand,1,lt) ;

I2 = normrnd(0,rand,1,lt) ;
E2 = normrnd(0,rand,1,lt) ;

EIc = normrnd(0,rand,1,lt) ; % common with cell and between cells

EI1 = normrnd(0,rand,1,lt) ; % common only within cell
EI2 = normrnd(0,rand,1,lt) ;


% waveforms in cell 1 and 2
Et1 = Ec+E1+EIc+EI1 ; 
Et2 = Ec+E2+EIc+EI2 ;

It1 = Ic+I1+EIc+EI1 ; 
It2 = Ic+I2+EIc+EI2 ; 

% what we can measure
measured_Et1 = sum(Et1.^2) ;
measured_It1 = sum(It1.^2) ;

measured_Et1Et2 = sum(Et1.*Et2) ; 
measured_It1It2 = sum(It1.*It2) ;

measured_Et1It2 = sum(Et1.*It2) ;
measured_Et2It1 = sum(Et2.*It1) ;

% what we want to solve for
actual_Et1It1(a) = sum(Et1.*It1) ; 
actual_Et1It1_Norm(a) = actual_Et1It1(a)/sqrt(measured_Et1*measured_It1) ;

% theoretical bounds
WeakUpperBound_Et1It1(a) = min(measured_Et1Et2,measured_It1It2) + min(measured_Et1-measured_Et1Et2,measured_It1-measured_It1It2) ;
UpperBound_Et1It1(a) = measured_Et1It2 + min(measured_Et1-measured_Et1Et2,measured_It1-measured_It1It2) ;
LowerBound_Et1It1(a) = measured_Et1It2 ;

WeakUpperBound_Et1It1_Norm(a) = WeakUpperBound_Et1It1(a)/sqrt(measured_Et1*measured_It1) ;
UpperBound_Et1It1_Norm(a) = UpperBound_Et1It1(a)/sqrt(measured_Et1*measured_It1) ;
LowerBound_Et1It1_Norm(a) = LowerBound_Et1It1(a)/sqrt(measured_Et1*measured_It1) ;

end

figure
plot(actual_Et1It1,actual_Et1It1,'k-')
hold on
plot(actual_Et1It1,WeakUpperBound_Et1It1,'b*')
plot(actual_Et1It1,UpperBound_Et1It1,'r*')
plot(actual_Et1It1,LowerBound_Et1It1,'g*')
xlabel('actual Et1It1')
ylabel('bounds')
legend('unity','weak upper','upper','lower')

figure
plot(actual_Et1It1_Norm,actual_Et1It1_Norm,'k-')
hold on
plot(actual_Et1It1_Norm,WeakUpperBound_Et1It1_Norm,'b*')
plot(actual_Et1It1_Norm,UpperBound_Et1It1_Norm,'r*')
plot(actual_Et1It1_Norm,LowerBound_Et1It1_Norm,'g*')
xlabel('actual Et1It1')
ylabel('bounds')
legend('unity','weak upper','upper','lower')

