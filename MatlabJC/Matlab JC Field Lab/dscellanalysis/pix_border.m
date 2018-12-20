function [mm,nn] = pix_border(m,n) 
% m: x_position
% n: y_postion
%Calculates 4 coordinate values for each pixel

    for i = 1:length(m)
        mm((i-1)*4+1) = m(i)+0.5;
        mm((i-1)*4+2) = m(i)+0.5;
        mm((i-1)*4+3) = m(i)-0.5;
        mm(i*4) = m(i)-0.5;
        nn((i-1)*4+1) = n(i)+0.5;
        nn((i-1)*4+2) = n(i)-0.5;
        nn((i-1)*4+3) = n(i)+0.5;
        nn(i*4) = n(i)-0.5;
    end
    mm = mm'; nn = nn';
end