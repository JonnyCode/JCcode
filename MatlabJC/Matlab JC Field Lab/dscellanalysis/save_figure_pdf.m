function[] = save_figure_pdf(pathname, filename, h)

%Function saves figure (h) as a pdf in folder (pathname) with given name (filename)

addpath(pathname)
fullpathname = [pathname filename];
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(h,'-dpdf', fullpathname);

end