% return x axis position corresponding to input value of cumulative gaussian specified
function XPosition = CumulativeGaussInv(Target, mu, sigma, scale)

XPosition = erfinv(2 * (Target / scale - 1/2)) * sigma * sqrt(2) + mu;

% check
% fprintf(1, '%d %d\n', normcdf(XPosition, mu, sigma), Target);
