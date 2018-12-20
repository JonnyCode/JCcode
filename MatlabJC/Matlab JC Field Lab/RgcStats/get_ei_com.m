function center_ei = get_ei_com(datarun, id, distance)

% XY 

frame = datarun.ei.nrPoints + datarun.ei.nlPoints + 1;
elec = size(datarun.ei.position, 1);
idx = get_cell_indices(datarun, id);
center_ei = zeros(length(id),2);
for i = 1:length(id)
    distance_temp = distance;
    full_elec_n = sum(1:distance_temp)*6+1; % Assume all arrays have hexagonal-arranged electrodes
    ei = datarun.ei.eis{idx(i)};
    ei = ei';
    [~,I] = max(abs(ei(:)));
    elec_n = ceil(I/frame);
    elecs_n = get_ei_neighbors(elec_n, elec, distance_temp);
    while(length(elecs_n) < full_elec_n)
        distance_temp = distance_temp - 1;
        elecs_n = get_ei_neighbors(elec_n, elec, distance_temp);
        full_elec_n = sum(1:distance_temp)*6+1;
    end
    points = datarun.ei.position(elecs_n,:);
    mass = max(abs(ei(:,elecs_n)));
    center_ei(i,:) = centroid(points, mass);
end

end
