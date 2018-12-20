%Function that normalizes by the value at direction = 0 degrees
%Inputs are rho (the values) and theta(the angles)
%Function returns the normalized radius and angle matrix after removing the rows that divide
%by zeros and have Nans and Inf

function[rho theta] = normalize_totalspikes(rho, theta)
for j = 1:length(rho)
    norm = [];
    norm = sum(rho{j,1}')';
    norm = repmat(norm, 1, size(rho{j,1},2));
    rho{j,1} = rho{j,1}./norm;
    [rho{j,1} theta{j,1}] = exciseRows(rho{j,1}, theta{j,1});

end
end