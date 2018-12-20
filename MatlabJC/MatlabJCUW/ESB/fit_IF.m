%seek to fit:
%C dV/dt = I_ion + I_app
%rewrite
%dV / dt = f(V) + I_app /C
%so we must find both f(V) and C!

load cell_data.mat

dt=tlist(2)-tlist(1);
dVdtlist = diff(Vlist)/dt;

%truncate all lists:  drop last point, so length matches dVdtlist
Vlist=Vlist(1:end-1);
gelist=gelist(1:end-1);
gilist=gilist(1:end-1);
tlist=tlist(1:end-1);
Iapplist=Iapplist(1:end-1);

%compute or load Iapplist
%Iapplist = gelist.*(Ve-Vlist) + gelist.*(Vi-Vlist) ;


%Infer or simply enter value of capacitance
%C=1   %actually is loaded

f_of_V_list = dVdtlist - Iapplist/C;  


%march through the list of spike times.  Label timesteps w/in certain
%ranges of every spike (moving forward) as to what window they fall into.

win1=1;
win2=1;
win3=1;
win4=1;
win5=Inf;

win_labels=zeros(1,length(tlist));

tpad=2*dt;

for j=1:length(tslist)-1
    j
    ts=tslist(j);
    tsnext=tslist(j+1);
    win_labels(find( ts+tpad < tlist & tlist < min(ts+win1,tsnext-tpad) ))=1;
    win_labels(find( ts+win1 < tlist & tlist < min(ts+win1+win2,tsnext-tpad) ))=2;
    win_labels(find( ts+win1+win2 < tlist & tlist < min(ts+win1+win2+win3,tsnext-tpad) ))=3;
    win_labels(find( ts+win1+win2+win3 < tlist & tlist < min(ts+win1+win2+win3+win4,tsnext-tpad) ))=4;
    win_labels(find( ts+win1+win2+win3+win4 < tlist & tlist < min(ts+win1+win2+win3+win4+win5,tsnext-tpad) ))=5;

end   
    
Vlist_win1=Vlist(find(win_labels==1));
Vlist_win2=Vlist(find(win_labels==2));
Vlist_win3=Vlist(find(win_labels==3));
Vlist_win4=Vlist(find(win_labels==4));
Vlist_win5=Vlist(find(win_labels==5));

f_of_V_list_win1=f_of_V_list(find(win_labels==1));
f_of_V_list_win2=f_of_V_list(find(win_labels==2));
f_of_V_list_win3=f_of_V_list(find(win_labels==3));
f_of_V_list_win4=f_of_V_list(find(win_labels==4));
f_of_V_list_win5=f_of_V_list(find(win_labels==5));



figure
subplot(511)
plot(Vlist_win1,f_of_V_list_win1,'.')
title('window 1')
hold on;
%superpose exact answer IN RED, to check:
plot(Vlist, 1/C*(-Vlist + Vlist.^2),'r.','LineWidth',4)
xlabel('V'); ylabel('f(V)')


subplot(512)
plot(Vlist_win2,f_of_V_list_win2,'.')
title('window 2')
hold on;
%superpose exact answer IN RED, to check:
plot(Vlist, 1/C*(-Vlist + Vlist.^2),'r.','LineWidth',4)
xlabel('V'); ylabel('f(V)')


subplot(513)
plot(Vlist_win3,f_of_V_list_win3,'.')
title('window 3')
hold on;
%superpose exact answer IN RED, to check:
plot(Vlist, 1/C*(-Vlist + Vlist.^2),'r.','LineWidth',4)
xlabel('V'); ylabel('f(V)')


subplot(514)
plot(Vlist_win4,f_of_V_list_win4,'.')
title('window 4')
hold on;
%superpose exact answer IN RED, to check:
plot(Vlist, 1/C*(-Vlist + Vlist.^2),'r.','LineWidth',4)
xlabel('V'); ylabel('f(V)')


subplot(515)
plot(Vlist_win5,f_of_V_list_win5,'.')
title('window 5')
hold on;
%superpose exact answer IN RED, to check:
plot(Vlist, 1/C*(-Vlist + Vlist.^2),'r.','LineWidth',4)
xlabel('V'); ylabel('f(V)')


figure
plot(Vlist,f_of_V_list,'.')
axis([-Inf Inf -5 5])
hold on;
%superpose exact answer IN RED, to check:
plot(Vlist, 1/C*(-Vlist + Vlist.^2),'r.','LineWidth',4)
xlabel('V'); ylabel('f(V)')





















