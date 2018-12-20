function[T R Unew Vnew] = polar_plots_one(rho,theta,U,V,num,cellnum)

%Plots polar plot for one cell for all speeds, also plots vector average or
%vector max
%Input rho and theta and preferred vector U and V of all cells for all stimuli for all speeds

axes_handle = [];
ylim = zeros(1, length(U));
T = cell(length(U),1);
R = cell(length(U),1);
Unew = zeros(length(U),1);
Vnew = zeros(length(U),1);
for i = 1:length(U)
    [thetanew, rhonew, ind] = deal([]);
    rhonew = rho{i,1}(cellnum,:);
    thetanew = theta{i,1}(cellnum,:);
    [thetanew, ind] = sort(thetanew);
    rhonew = rhonew(ind);
    Unew(i,1) = U{i,1}(1,cellnum);
    Vnew(i,1) = V{i,1}(1,cellnum);
    T{i,1} = thetanew;
    R{i,1} = rhonew;
    ax(i) = subplot(1,length(U),i); %Plots different speeds on different subplot
    axes_handle = [axes_handle ax(i)];
    polar(thetanew, rhonew); %Polar plot
    ylim(1,i) = max(rhonew);
    hold on;
    h = compass(Unew(i,1),Vnew(i,1), 'r'); %Vector average or vector sum plot
    set(h,'linewidth',3) 
    titlechar = ['Polar and Vector Average Plot for cell ' num2str(cellnum) ' of temporal period: ' num2str(num(i))];
    title(titlechar, 'Color', 'k', 'FontWeight', 'bold' , 'FontSize', 18 ,'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Right');
    hold off
end

linkaxes(axes_handle,'xy'); %Link axes to have same x and y limits
[C, I] = max(ylim);
set(ax(I),'xlimmode','auto');
set(ax(I),'ylimmode','auto');

end



