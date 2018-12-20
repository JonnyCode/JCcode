% this simulates the error in the measured noise correlations caused during
% simultaneous coductance recordings when the holding voltage is off.  This
% is equivilent to equation 2.7 in "Correlated noise..." paper.

c=[0:.1:1] ; % correlation
v2=[-50.0005:10:0] ; % voltage hold

Me_var = 10 ; % variance exc g (var(Ge)+var(s))
Mi_var = 10 ;

Ee=0 ; % revesal exc
Ei=-60;

v1=Ei; %

t=100000 ;

s_orig = normrnd(10,1,1,t) ; % conductances
Ge_orig = normrnd(10,1,1,t) ;
Gi_orig = normrnd(10,1,1,t) ;


for trial = 1:2 ;
   

for c_round = 1:length(c) ; % for each new correlation value

    
    for v_round = 1:length(v2) ; % for each new hold     
    
    S_var = c(c_round)*sqrt(Me_var*Mi_var) ; % variance of shared g
    
    s = s_orig*sqrt(S_var) ; % g actual
    Ge = Ge_orig*sqrt(Me_var-S_var) ;
    Gi = Gi_orig*sqrt(Mi_var-S_var) ;

    Me = ((Ge+s)*(v1-Ee)+(Gi+s)*(v1-Ei))/(Ei-Ee) ; % g measured
    Mi = ((Ge+s)*(v2(v_round)-Ee)+(Gi+s)*(v2(v_round)-Ei))/(Ee-Ei) ;
    
    if trial == 1 ;
       Mi = ((Ge+s)*(v2(v_round)-Ee)+(0)*(v2(v_round)-Ei))/(Ee-Ei) ; % g with inh blocked
    end
    
    AcMi(c_round,v_round) = max(xcov(Mi)) ; % variance of measured g
    AcMe(c_round,v_round) = max(xcov(Me)) ;
    
    Corr_actual = xcov(Ge+s,Gi+s) ; % actual g correlation numerator
    cna(c_round,v_round) = Corr_actual(find(abs(Corr_actual)==max(abs(Corr_actual)))) ;
    
    Corr_measured = xcov(Me,Mi) ; % measured g correatlion numerator
    cnm(c_round,v_round) = Corr_measured(find(abs(Corr_measured)==max(abs(Corr_measured)))) ;
    
    cda(c_round,v_round) = max(sqrt(xcov(Ge+s).*xcov(Gi+s))) ; % denominator actual
    cdm(c_round,v_round) = max(sqrt(xcov(Me).*xcov(Mi))) ; % denominator measured
    
    ca(c_round,v_round) = cna(c_round,v_round)/cda(c_round,v_round) ; % actuall corr coef
    cm(c_round,v_round) = cnm(c_round,v_round)/cdm(c_round,v_round) ; % measured corr coef
    
    
    end
end

if trial == 1 ;
    cnmInh = cnm ;
    AcMiInh = AcMi ;
    
    clearvars -except cnmInh AcMiInh c v1 v2 Me_var Mi_var Ee Ei t s_orig Ge_orig Gi_orig;    
else    
%     cnmCorrected = cnm - cnmInh ; % correcting the wrong way
%     cmCorrected = cnmCorrected./cdm ;
    
    d = sqrt(AcMiInh./AcMe) ; % corrected the right way
   
    % % correction for linear i-v plots
    %cnCorrected = (cnm - cnmInh)./(1-d) ; % numerator corrected
    %cdCorrected = sqrt(((AcMi - (d.^2).*AcMe + 2*d.*(cnm - cnmInh))./(1-d).^2).*AcMe) ; % denominator corrected
    
    % corrections for nonlinear i-v plots
    cnCorrected = (cnm - cnmInh) ; % numerator corrected
    cdCorrected = sqrt((AcMi - (d.^2).*AcMe + 2*d.*(cnm - cnmInh)).*AcMe) ; % denominator corrected
    
    cCorrected = cnCorrected./cdCorrected ; % corrected correlation
    
    
    figure
    plot(cnm,cna,'*-')
    hold on
    plot(cnm,cnCorrected,'o--')
    
    figure
    plot(cdm,cda,'-*')
    hold on
    plot(cdm,cdCorrected,'o--')
    
    figure
    plot(cm,ca,'*-')
    hold on
    plot(cm,cCorrected,'o--')    
    ylabel('actual correlation coef')
    xlabel('measured correlation coef')
    
%     legend('Ee=-20', 'Ee=-10','Ee=0','Ee=-20,with corrections','Ee=-10,with corrections','Ee=0,with corrections')
end

end

figure
plot(cm,ca,'*-')
hold on
plot(cm,ca,'o--')
ylabel('actual correlation coef')
xlabel('measured correlation coef')


figure
plot((v2-Ee)/60,v2-Ee)
xlabel('fraction of variance Exc reduced at assumed Ee as compared to Ei')
ylabel('voltage driving force equivelent')



% figure
% plot(c,error,'*-')
% 
figure
plot(v2-Ee,errorcn)
% 
figure
plot(v2-Ee,errorcd)
% 
figure
plot(v2-Ee,error)
% 
% figure
% plot(v2-Ee,cm(6,:))
% 
% figure
% plot(errorcn,errorcd,'*')
% 
% figure
% plot(errorcn(6,:),errorcd(6,:),'*')
% 
% figure
% plot(errorcn(:,end-1),errorcd(:,end-1),'*-')
% hold on
% plot(errorcn(1,end-1),errorcd(1,end-1),'o')

% figure 
% plot(cnm,v2-Ee,'-*')
% 
% figure
% plot(ca(:,end-1),cm(:,end-1),'b-*')
% hold on
% plot(ca(:,end-1),ca(:,end-1),'k')



