%Calculation of Parameters from Autocorrelation Function
function [ACFParams] = acf_parameters(autocorrel, bins)

% Input: Auto Correlation matrix (rows - Acf values, columns - each cell) 
% Returns structure with parameters calculated from each cell's autocorrelation function in this order:

% 1. Standard Deviation
% 2. Variance
% 3. Mean
% 4. Minimum Value
% 5. Minimum Value Time
% 6. Peak Value
% 7. Peak Value Time
% 8. Norm
% 9. Absolute Area Under Curve
%10. Total Area Under Curve
%11. 

%26. Minimum Gradient
%27. Minimum Gradient Time
%28. Peak Gradient
%29. Peak Gradient Time
%30. Minimum Gradient Time*Peak Gradient Time
%31. Minimum Gradient Time+Peak Gradient Time
%32. Peak Gradient Time - Minimum Gradient Time
%33. Peak Gradient Time / Minimum Gradient Time
%34. Minimum Gradient*Peak Gradient 
%35. Minimum Gradient+Peak Gradient 
%36. Peak Gradient - Minimum Gradient
%37. Peak Gradient / Minimum Gradient
%38. Abs(Minimum Gradient+Peak Gradient)









ACFParams = struct;

ACFParams.stddev = std(autocorrel);
ACFParams.var = var(autocorrel);
ACFParams.ave = mean(autocorrel);
[ACFParams.minval mintim] =  min(autocorrel); %Minimum value in ACF doesn't really mean anything, can be anywhere (before or after peak)
ACFParams.mintim = bins(1,mintim);
[ACFParams.maxval maxtim] =  max(autocorrel);
ACFParams.maxtim = bins(1,maxtim);

for i = 1:size(autocorrel, 2)
  ACFParams.norm(1, i) = norm( autocorrel(:,i)); %Calculate norm (magnitude) for all acf
  ACFParams.areaabs(1,i) = trapz(abs(bins), abs(autocorrel(:,i))); %Calculate Area Under Curve for all acf
  ACFParams.area(1,i) = trapz(abs(bins), autocorrel(:,i));
end 

[ACFParams.mingradval ACFParams.mingradtim] =  min(diff(autocorrel)); 
[ACFParams.maxgradval ACFParams.maxgradtim] =  max(diff(autocorrel)); 

ACFParams.mingradtmaxgradtmult = ACFParams.mingradtim.*ACFParams.maxgradtim;
ACFParams.mingradtmaxgradtadd = ACFParams.mingradtim+ACFParams.maxgradtim;
ACFParams.mingradtmaxgradtsub = ACFParams.maxgradtim - ACFParams.mingradtim;
ACFParams.mingradtmaxgradtdiv = ACFParams.maxgradtim./ACFParams.mingradtim;

ACFParams.mingradmaxgradmult = ACFParams.mingradval.*ACFParams.maxgradval;
ACFParams.mingradmaxgradadd = ACFParams.mingradval+ACFParams.maxgradval;
ACFParams.mingradmaxgradsub = ACFParams.maxgradval - ACFParams.mingradval;
ACFParams.mingradmaxgraddiv = ACFParams.maxgradval./ACFParams.mingradval;
ACFParams.mingradmaxgradabsadd = abs(ACFParams.mingradval) + abs(ACFParams.maxgradval);




    
end