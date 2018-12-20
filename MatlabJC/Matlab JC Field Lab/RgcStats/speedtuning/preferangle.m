function[prefer] = preferangle(rho, angle, theta)

%Function returns spike number of direction with least firing

% Sneha Ravi 
% Last revision: 12-18-2012

prefer = cell(size(rho));%
for i = 1:size(rho,1)
    for j = 1:size(rho,2)
        rho1 = rho{i,j};
        theta1 = theta{i,j};
        angle1 = angle{i,j};
        angle1 = mod(angle1, 2*pi);
        angle1 = repmat(angle1', 1, size(theta1, 2));
        [~, preferi] = min(abs(angle1-theta1), [], 2);
        for cc = 1:size(theta1, 1)
            prefer{i,j}(cc) = rho1(cc, preferi(cc));
        end
    end
end