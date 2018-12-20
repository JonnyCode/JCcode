clear all
close all

M_ee = 1.25
M_ii = -1
M_ei = -1
M_ie = 1
gamma_e = -10
gamma_i = 10
tau_e = 10

V_e = 40
V_i = 30
Ve_trace=[V_e]
Vi_trace=[V_i]

tau_i = 80

for round= 1:1000

V_e = V_e + (1/tau_e)*(-V_e +(max(((M_ee*V_e))+((M_ei*V_i))-(gamma_e),0))) ;

Ve_trace=[Ve_trace V_e] ;

V_i = V_i + (1/tau_i)*(-V_i + (max(((M_ii*V_i))+((M_ie*V_e))-(gamma_i),0))) ;

Vi_trace=[Vi_trace V_i] ;

end

figure, plot(Vi_trace)
xlabel('time')
ylabel('V_i Hz')

figure, plot(Ve_trace)
xlabel('time')
ylabel('V_e Hz')

figure, plot(Ve_trace,Vi_trace)
xlabel('V_e')
ylabel('V_i')
