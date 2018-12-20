figure

plot(ForIgor.SurroundFactHistDb29(1,:),ForIgor.SurroundFactHistDb29(2,:)/(ForIgor.SurroundFactHistDb29(2,end)),'k')
hold on
plot(ForIgor.SurroundFactHistDb30(1,:),ForIgor.SurroundFactHistDb30(2,:)/(ForIgor.SurroundFactHistDb30(2,end)),'k')
plot(ForIgor.SurroundFactHistDb31(1,:),ForIgor.SurroundFactHistDb31(2,:)/(ForIgor.SurroundFactHistDb31(2,end)),'k')
plot(ForIgor.SurroundFactHistDb32(1,:),ForIgor.SurroundFactHistDb32(2,:)/(ForIgor.SurroundFactHistDb32(2,end)),'r')
plot(ForIgor.SurroundFactHistDb33(1,:),ForIgor.SurroundFactHistDb33(2,:)/(ForIgor.SurroundFactHistDb33(2,end)),'r')
plot(ForIgor.SurroundFactHistDb34(1,:),ForIgor.SurroundFactHistDb34(2,:)/(ForIgor.SurroundFactHistDb34(2,end)),'r')
legend('control','','','+plp')
xlabel('surround factor')
ylabel('fraction of cells')


figure
plot(ForIgor.ZeroTimeHistDb29(1,:),ForIgor.ZeroTimeHistDb29(2,:)/(ForIgor.ZeroTimeHistDb29(2,end)),'k')
hold on
plot(ForIgor.ZeroTimeHistDb30(1,:),ForIgor.ZeroTimeHistDb30(2,:)/(ForIgor.ZeroTimeHistDb30(2,end)),'k')
plot(ForIgor.ZeroTimeHistDb31(1,:),ForIgor.ZeroTimeHistDb31(2,:)/(ForIgor.ZeroTimeHistDb31(2,end)),'k')
plot(ForIgor.ZeroTimeHistDb32(1,:),ForIgor.ZeroTimeHistDb32(2,:)/(ForIgor.ZeroTimeHistDb32(2,end)),'r')
plot(ForIgor.ZeroTimeHistDb33(1,:),ForIgor.ZeroTimeHistDb33(2,:)/(ForIgor.ZeroTimeHistDb33(2,end)),'r')
plot(ForIgor.ZeroTimeHistDb34(1,:),ForIgor.ZeroTimeHistDb34(2,:)/(ForIgor.ZeroTimeHistDb34(2,end)),'r')
legend('control','','','+plp')
xlabel('zero cross time (s)')
ylabel('fraction of cells')

figure
plot(ForIgor.DotHistDb29(1,:),ForIgor.DotHistDb29(2,:)/(ForIgor.DotHistDb29(2,end)),'k')
hold on
plot(ForIgor.DotHistDb30(1,:),ForIgor.DotHistDb30(2,:)/(ForIgor.DotHistDb30(2,end)),'k')
plot(ForIgor.DotHistDb31(1,:),ForIgor.DotHistDb31(2,:)/(ForIgor.DotHistDb31(2,end)),'k')
plot(ForIgor.DotHistDb32(1,:),ForIgor.DotHistDb32(2,:)/(ForIgor.DotHistDb32(2,end)),'r')
plot(ForIgor.DotHistDb33(1,:),ForIgor.DotHistDb33(2,:)/(ForIgor.DotHistDb33(2,end)),'r')
plot(ForIgor.DotHistDb34(1,:),ForIgor.DotHistDb34(2,:)/(ForIgor.DotHistDb34(2,end)),'r')
legend('control','','','+plp')
xlabel('degree of transients')
ylabel('fraction of cells')
