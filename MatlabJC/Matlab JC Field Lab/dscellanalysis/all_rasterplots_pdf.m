function[] = all_rasterplots_pdf(NumSpikesCell,StimComb,datarun, rho,theta,U,V,num, chos, save,plots)

%cellnum = DSCELLS(1,:);

for cellnum = 1: length(chos(1,:))
% addpath('/Users/sravi/matlab/DS cell analysis/2012-10-15-0/Data002PolarRasterPlots');
% pathname = ['/Users/sravi/matlab/DS cell analysis/2012-10-15-0/Data002PolarRasterPlots/'  num2str(chos(1,cellnum)) '-' num2str(datarun.cell_ids(chos(1,cellnum)))];
[T R Unew Vnew] = polar_plots_one(rho,theta,U,V,num,chos(1,cellnum));
close;
if(plots == 1)
    raster_plot(StimComb,datarun,num,T,R,Unew, Vnew, chos(1,cellnum));
elseif(plots == 2)
    raster_plots(StimComb,datarun,num,T,R,Unew, Vnew, chos(1,cellnum));
end
if(save)
    save_figure_pdf('/Analysis/sravi/Rat/WildType/2012-10-31-0/DSPlots/', [num2str(chos(1,cellnum)) '-' num2str(datarun.cell_ids(chos(1,cellnum)))], gcf);
end

end

end


