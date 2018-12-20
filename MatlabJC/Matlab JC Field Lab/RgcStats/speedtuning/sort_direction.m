function ds_out = sort_direction(ds_in)

% ds_in: structure contain all ds information
% xyao 
% 2014-08-12

for j = 1:size(ds_in.rho, 1)
    for i = 1:size(ds_in.rho, 2)
        for k = 1:size(ds_in.rho, 3)
            [theta_seq, I_seq] = sort(ds_in.theta{j, i, k}(1, :));
            r = ds_in.rho{j, i, k};
            R = ds_in.RHO{j, i, k};
            rho_seq = r(:, I_seq);
            RHO_seq = R(:, I_seq);
            ds_in.rho{j, i, k} = rho_seq;
            ds_in.RHO{j, i, k} = RHO_seq;
            ds_in.theta{j, i, k} = repmat(theta_seq, size(ds_in.theta{j, i, k}, 1), 1);
        end
    end
ds_out = ds_in;
end