function[null] = nullopp(rho, angle, theta)

%Function returns spike number of direction nearest to that opposite direction of vector ave/sum/max

null = cell(length(rho), 1);

for i = 1:length(rho)
    b = angle{i,1}; 
    b(b<0) = b(b<0)+(2*pi); %b is the angles from 0 to 360 degrees
    c = b <= pi;
    d = b > pi;
    b(c) = b(c)+ pi;
    b(d) = b(d)- pi;
    thetnew = sort(theta{i,1}(1,:));
    bnew = [];
    thetnew(1,length(thetnew)+1) = 2*pi;
    bnew = interp1(thetnew,thetnew,b,'nearest'); %null angle
    bnew(bnew == max(thetnew)) = min(thetnew);
    bnew = bnew';
    bnew = repmat(bnew, 1, size(theta{i,1},2));
    yo = (bnew==theta{i,1});
    rhonew = rho{i,1}';
    yo = yo';
    null{i,1} = rhonew(yo)';
end

end
    