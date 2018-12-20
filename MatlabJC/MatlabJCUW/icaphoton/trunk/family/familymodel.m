function y = familymodel(k,I)
% y = familymodel(k,I) 			
% function that is used to fit the Sensitivity vs. log of flash
% intensity.
% I is the flash intensity.
y = 1 - exp(-k*I);
