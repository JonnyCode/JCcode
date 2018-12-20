cellids = datarun000.cell_ids;

for a = 1:length(cellids)
    off_cell_sta = get_sta(datarun000, cellids(1,a));
    off_cell_sta = squeeze(off_cell_sta);
    sta_dims = size(off_cell_sta);
    off_cell_sta = reshape(off_cell_sta, prod(sta_dims(1:2)), 30);
    % perform singular value decomposition
    [space_vecs, singular_values, time_vecs] = svd(off_cell_sta, 0);
    sep_index(1,a) = singular_values(1) ./ sum(diag(singular_values));
%     spatial_singular_vector = reshape(space_vecs(:,1), sta_dims(1:2));
%     figure(4)
%     imagesc(spatial_singular_vector); colormap gray;
%     title('first spatial singular vector: inseparable')
end

plot(1:length(cellids), sep_index, 'o');


cellids = datarun000.cell_ids;
%cellids = Cs;
for a = 1:length(cellids)
    off_cell_sta = get_sta(datarun000, cellids(1,a));
    off_cell_sta = squeeze(off_cell_sta);
    off_cell_stax = sum(off_cell_sta,1);
    sta_dimsx = size(off_cell_stax);
    off_cell_stax = reshape(off_cell_stax, prod(sta_dimsx(1:2)), 30)';
    
    off_cell_stay = sum(off_cell_sta,2);
    sta_dimsy = size(off_cell_stay);
    off_cell_stay = reshape(off_cell_stay, prod(sta_dimsy(1:2)), 30)';
    figure(1);
    contour(off_cell_stax);
    title([num2str(cellids(1,a)), '-', num2str(sep_index1(1,a))]);
     figure(2);
   contour(off_cell_stay);
    pause;
    close all;
end



cellids = datarun000.cell_ids;

for a = 1:length(cellids1)
    off_cell_sta = get_sta(datarun000, cellids1(1,a));
    off_cell_sta = squeeze(off_cell_sta);
    off_cell_sta = sum(off_cell_sta,2);
    sta_dims = size(off_cell_sta);
    off_cell_sta = reshape(off_cell_sta, prod(sta_dims(1:2)), 30)';
    plot(off_cell_sta);
%     contour(off_cell_sta);figure(gcf);
    title([num2str(cellids1(1,a)), '-', num2str(sep_index1(1,a))]);
    pause;
end

ds = [1595 2042 301  3736 4069 438  4846 5632 6321 708 1683 2449 333  3815 4157 467  4985 5702 6332 7308 1685 257  3512 3842 4173 4731 5148 5719 6751 7475 1097 1895 2898 3636 3934 4353 4789 5150 6229 6797 766];
dsind = get_cell_indices(datarun000,ds);
dsind = sort(dsind);
sep_index(1,dsind)


Cs = [257 301 333 424 438 573 708 857 1084 1097 1114 1130 1191 1352 1546 1564 1595 1683 1685 1863 2042 2132 2206 2236 2299 2311 2449 2795 2898 3046 3079 3364 3512 3586 3636 3648 3662 3679 3736 3767 3815 3841 3842 3859 3889 3934 4022 4025 4069 4157 4173 4188 4279 4294 4353 4368 4456 4622 4731 4788 4789 4834 4846 4940 4955 4985 5042 5043 5131 5134 5150 5632 5645 5719 6049 6064 6181 6229 6321 6332 6366 6726 6751 6797 6935 7040 7052 7054 7082 7278 7308 7310 7431 7475 7610];

intersect(ds, Cs)



cellids = datarun000.cell_ids;
cellids = sort(Cs);
for a = 1:length(cellids)
    rf = get_rf(datarun000,cellids(1,a));
    off_cell_sta = get_sta(datarun000, cellids(1,a));
    off_cell_sta = squeeze(off_cell_sta);
    ctr = round(rf_center(datarun000,cellids(1,a)));
    rad = 10;
    xrng = max(ctr(1)-rad,1):min(ctr(1)+rad,size(off_cell_sta,2)); %80
    yrng = max(ctr(2)-rad,1):min(ctr(2)+rad,size(off_cell_sta,1)); %40
    off_cell_sta = off_cell_sta(yrng, xrng, :);
    off_cell_stax = sum(off_cell_sta,1);
    sta_dimsx = size(off_cell_stax);
    off_cell_stax = reshape(off_cell_stax, prod(sta_dimsx(1:2)), 30)';
    
    off_cell_stay = sum(off_cell_sta,2);
    sta_dimsy = size(off_cell_stay);
    off_cell_stay = reshape(off_cell_stay, prod(sta_dimsy(1:2)), 30)';
    figure(1);
    contourf(off_cell_stax);
    %plot(mean(off_cell_stax));
    title([num2str(cellids(1,a)), '-', num2str(sep_index1(1,a))]);
%     figure(2);
%     contour(off_cell_stay);
    pause;
    close all;
end


for a = 1:length(cellids)
    for b = 1:30
         com = rf_com(datarun000.stas.stas{cellids(1,a),1}(:,:,:,b));
         if(isempty(com))
             comx(a,b) = 0;
             comy(a,b) = 0;
         elseif (isnan(com))
             comx(a,b) = 0;
             comy(a,b) = 0;
         else
             comx(a,b) = com(1,1);
             comy(a,b) = com(1,2);
         end
    end
end

for c = 1:length(cellids)
    xcoor = [];
    ycoor = [];
     e = 1;
     f = 1;
    for d = 1:size(comx,2)
        if (comx(c,d) ~= 0)
            xcoor(1,e) = comx(c,d);
            e = e+1;       
        end
        if (comy(c,d) ~= 0)
            ycoor(1,f) = comy(c,d);
            f = f+1;
        end
    end
    figure();
    line(xcoor, ycoor);
    axis([0 80 0 40])
    title(num2str(ds(1,c)));
    pause;
    close;
end
    
[x,y] = meshgrid(comx(c,:), comy(c,:));



