% saving figures from matlab for ppt presenations

print(gcf,'-dpng','-r300',['ORNfig1d2', '.png']); % saves a high resolution image
set(gca,'LineWidth',1) ; % sets thick axes
set(gca,'box','off') ; % turns off box
set(gca,'TickLength',[.02,.02]) ; % tick length (normalized units [2d, 3d]

set('LineWidth',2)