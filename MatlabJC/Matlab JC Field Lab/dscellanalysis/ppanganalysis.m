  plot(mag{1,1}, mag{2,1}, '+')
    plot(mag{1,1}, mag{3,1}, '+')
    plot(mag{3,1}, mag{2,1}, '+')
     plot(mag{2,1}, mag{3,1}, '+')

       plot(magin{1,1}, magin{2,1}, '+')

  ginput(1);

Y = get(get(gca,'Children'),'YData');

X = get(get(gca,'Children'),'XData');

h = impoly;

accepted_pos = wait(h);

hold on;

drawPolygon(accepted_pos(:,1), accepted_pos(:,2));

in = inpolygon(X, Y, accepted_pos(:,1), accepted_pos(:,2));

plot(X(in),Y(in),'r+');

cid = 1:length(X);

chos = cid(in);



UU{1,1} = U{1,1}(1, chos(1,:))
UU{2,1} = U{2,1}(1, chos(1,:))
UU{3,1} = U{3,1}(1, chos(1,:))

VV{1,1} = V{1,1}(1, chos(1,:))
VV{2,1} = V{2,1}(1, chos(1,:))
VV{3,1} = V{3,1}(1, chos(1,:))

magmag{1,1} = mag1{1,1}(1, chos(1,:))
magmag{2,1} = mag1{2,1}(1, chos(1,:))
magmag{3,1} = mag1{3,1}(1, chos(1,:))

figure(1) 
polar_plots_all(UU,VV,num, magmag);
 
 plot(angle{2,1}(1, chos(1,:)), 1, '+')
 
   plot(anglein{1,1}(1, chos(1,:)), anglein{2,1}(1, chos(1,:)), '+')
     ginput(1);

Y = get(get(gca,'Children'),'YData');

X = get(get(gca,'Children'),'XData');

h = impoly;

accepted_pos = wait(h);

hold on;

drawPolygon(accepted_pos(:,1), accepted_pos(:,2));

in = inpolygon(X, Y, accepted_pos(:,1), accepted_pos(:,2));

plot(X(in),Y(in),'r+');

cidcid = 1:length(X);

choschos = cid(in);
   
   
   

 figure(2)
  plot(angle{1,1}(1, chos(1,:)), angle{2,1}(1, chos(1,:)), '+')
   figure(3)
  plot(angle{1,1}(1, chos(1,:)), angle{3,1}(1, chos(1,:)), '+')
   figure(4)
  plot(angle{2,1}(1, chos(1,:)), angle{3,1}(1, chos(1,:)), '+')
  
    plot(angle{2,1}(1, chos(1,:)), angle{3,1}(1, chos(1,:)), '+')
        plot(mag{2,1}(1, chos(1,:)), mag{3,1}(1, chos(1,:)), '+')
                plot(dsindex{2,1}(chos(1,:),1), dsindex{3,1}(chos(1,:),1), '+')

 
%data04 or 05%%%%%%%%%%%%%%%                
UU{1,1} = U{1,1}(1, chos(1,:))
UU{2,1} = U{2,1}(1, chos(1,:))
VV{1,1} = V{1,1}(1, chos(1,:))
VV{2,1} = V{2,1}(1, chos(1,:))
magmag{1,1} = mag1{1,1}(1, chos(1,:))
magmag{2,1} = mag1{2,1}(1, chos(1,:))
polar_plots_all(UU,VV,num, magmag);
 
plot(angle{1,1}(1, chos(1,:)), angle{2,1}(1, chos(1,:)), '+')
plot(dsindex{1,1}(chos(1,:),1), dsindex{2,1}(chos(1,:),1), '+')
plot(mag{1,1}(1, chos(1,:)), mag{2,1}(1, chos(1,:)), '+')                
plot(angle{1,1}(1, chos(1,:)), angle{2,1}(1, chos(1,:)), '+')

   plot(angle{2,1}(1, chos(1,:)), angle{3,1}(1, chos(1,:)), '+')   
   plot(angle{1,1}(1, chos(1,:)), angle{2,1}(1, chos(1,:)), '+')

      plot(angle{2,1}(1, chos(1,:)), angle{3,1}(1, chos(1,:)), '+')   

plot(anglein{1,1}(1, chos(1,:)), anglein{2,1}(1, chos(1,:)), '+')
ginput(1);
Y = get(get(gca,'Children'),'YData');
X = get(get(gca,'Children'),'XData');
h = impoly;
accepted_pos = wait(h);
hold on;
drawPolygon(accepted_pos(:,1), accepted_pos(:,2));
in = inpolygon(X, Y, accepted_pos(:,1), accepted_pos(:,2));
plot(X(in),Y(in),'r+');
chos2 = chos(in);
chos2 = cellids(1,chos2)

chos2 = [172 184 247 364 402 427 480 597];

cellindic = chos2;
cellids = datarun.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun, datarun01, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun000, [257        1097        1683        1895        2898        3512        3736        4069        4157        4846       5632        6321        6751        6797])
plot_rf_summaries(datarun000, [333         438         467        1595        1685        2042        3636        3842        4353        4731        4985        5150        5702        7475])
plot_rf_summaries(datarun000,ds(1,F2F1>1))
figure(2)
plot_rf_portraits(datarun01,A)
plot_rf_portraits(datarun000,ds(1,F2F1>1))
figure(3)
plot_rfs(datarun000, c)
plot_rfs(datarun000,ds(1,F2F1>1))

figure(4)
plot_time_courses(datarun000, c, 1, 0)


cellindic = chos2;
cellids = datarun.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
A(2,:) = datarun04.cell_ids(1,(emptyCells==0));
figure(1)
plot_rf_summaries(datarun000, A(1,:))
figure(2)
plot_rf_portraits(datarun000,A(1,:))
figure(3)
plot_rfs(datarun000, A(1,:))
figure(4)
plot_time_courses(datarun000, A(1,:), 1, 0)
%%%%%%%%%%%%%%%%%

%%data02

UU{1,1} = U{1,1}(1, DSCELLS(1,:));
UU{2,1} = U{2,1}(1,  DSCELLS(1,:));
VV{1,1} = V{1,1}(1,  DSCELLS(1,:));
VV{2,1} = V{2,1}(1,  DSCELLS(1,:));
magmag{1,1} = mag1{1,1}(1,  DSCELLS(1,:));
magmag{2,1} = mag1{2,1}(1,  DSCELLS(1,:));

 polar_plots_all(UU,VV,num, magmag);
 figure(2)
 plot(angle1{1,1}(1, DSCELLS(1,:)), angle1{2,1}(1, DSCELLS(1,:)), '+')
 
 %%%%%%%%%%%%%%%%%%

