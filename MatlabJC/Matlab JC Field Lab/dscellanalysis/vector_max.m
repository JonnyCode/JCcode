function[U V angle mag] = vector_max(rho, theta)

%Function that identifies maximum spike number over all directions     sravi 12-18-2012

%Input: X and Y magnitudes of all cells for all temporal periods in all directions

%Output: Cartesian (U,V) and Polar (angle, mag) forms of calculated maximum vector


U = cell(length(rho), 1);
V = cell(length(rho), 1);
angle = cell(length(rho), 1); 
mag = cell(length(rho), 1); 
for i = 1:length(rho)
    rhonew = rho{i,1}';
    [rhonew ind] = max(rhonew);
    ind(2,1:length(ind)) = 1:length(ind);
    thetanew = theta{i,1}(ind(2,:), ind(1,:));
    thetanew = thetanew(1,:);
    [U{i,1} V{i,1}] = pol2cart(thetanew,rhonew);
    [angle{i,1},mag{i,1}] = cart2pol(U{i,1}, V{i,1});
end
end