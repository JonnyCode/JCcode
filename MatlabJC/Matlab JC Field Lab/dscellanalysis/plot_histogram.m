function[] = plot_histogram(StimComb, NumSpikesCell)



scatter(1:868, max(mag))

plot(num(1), mag{1,1}, num(2), mag{2,1})

[C I] = max(mag); %Calculate maximum vector average magnitude and the speed it corresponds to


hist(mag(1,1:868));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w');
hold on
hist(mag(2,1:868));

%Histogram / Plots of maximum vector average magnitude, still gives
%unimodal distribution
plot(I, C, 'o');
hist(C(1,1:868));figure(gcf);
hist(C(1,1:868), 60);
hist(C(1,1:868), 100);
hist(C(1,1:868), 200);


%Histogram plots of spike rates at each speed with varying bins, still gives
%unimodal distribution
hist(vecAve(1,1:868));figure(gcf);
hist(vecAve(1,1:868), 60)
hist(vecAve(1,1:868), 120)

hist(vecAve(2,1:868));figure(gcf);
hist(vecAve(2,1:868), 100);

%Plot values from one speed vs other speed, can differentiate better
%between DS and non-DS cells
plot(mag(1,:), mag(2,:), 'o')
hold on
plot((1:60), (1:60)) %straight linear line look above and below
hold off
plot(log(mag(1,:)), log(mag(2,:)), 'o')


%dot product

M = U{1,1}.*U{2,1};
N = V{1,1}.*V{2,1};
dotP = M+N;
scatter(1:868, dotP);
hist(dotP(1,1:868) ,15);
plot(dotP(1,1:868), '+');







new = log(dotP);

plot(real(new), 'o'); % Complex numbers are produced if dotP is negative
plot(imag(new), 'o'); %positive numbers have 0 Im part, negative numbers have Im part
hist(real(new))
% If A and B are perpendicular (at 90 degrees to each other), the result of the dot product will be zero, because cos(?) will be zero.
% If the angle between A and B are less than 90 degrees, the dot product will be positive (greater than zero), as cos(?) will be positive, and the vector lengths are always positive values.
% If the angle between A and B are greater than 90 degrees, the dot product will be negative (less than zero), as cos(?) will be negative, and the vector lengths are always positive values.

%a = find(dotP<=0)
a = find(dotP>0);
hist(dotP(a), 20)








%%%%%%%%%%%DS CELLS ANGLE AND SPEED CLASSIFICATION %%%%%%%%%%%
plot(mag{1,1},mag{2,1}, 'o');
plot(mag{1,1}(1,chos), mag{2,1}(1,chos), '+')
plot(dsindex{1,1},dsindex{2,1}, '+');
plot(dsindex{1,1}(chos,1), dsindex{2,1}(chos,1), '+')
plot(angle{1,1}, angle{2,1}, '+')
plot(angle{1,1}(1, chos), angle{2,1}(1, chos), '+')
hist(mag{1,1}(1,chos), 50)
hist(mag{2,1}(1,chos), 50)
hist(magmax(1,chos), 50)
hist(magave(1,chos), 50)
hist(angle{1,1}(1,chos), 50)
hist(angle{2,1}(1,chos), 50)
hist(dsindex{1,1}(chos,1), 50)
hist(dsindex{2,1}(chos,1), 50)



plot(angle{1,1},angle{2,1}, '+');

plot(mag{1,1},mag{2,1}, '+');

plot(dsindex{1,1},dsindex{2,1}, '+');

plot(1:868, dsindex{1,1}, '+') 

plot(mag{1,1},mag{2,1}, '+');

plot(dsindex{1,1}(chos,1), dsindex{2,1}(chos,1), '+')


plot(mag{1,1},mag{2,1}, '+');

ginput(1);

Y = get(get(gca,'Children'),'YData');

X = get(get(gca,'Children'),'XData');

h = impoly;

accepted_pos = wait(h);

hold on;

drawPolygon(accepted_pos(:,1), accepted_pos(:,2));

in = inpolygon(X, Y, accepted_pos(:,1), accepted_pos(:,2));

plot(X(in),Y(in),'r+');

cid = 1:length(X);

chos = cid(in);

ismember(chos, DSCELLS(1,:))
sum(ans)
ismember(chos, DSCELLS(2,:))
sum(ans)


chos = DSCELLS(1,:);
magDS = [];
angleDS = [];
%magDS = cell(2,1);
magDS(1,:) = mag{1,1}(1,chos); %32 or fast
magDS(2,:) = mag{2,1}(1,chos); %256 or slow
angleDS(1,:) = angle{1,1}(1,chos);
angleDS(2,:) = angle{2,1}(1,chos);

angleDS(1,(angleDS(1,:) <0)) = angleDS(1,(angleDS(1,:) <0)) + 2*pi;
angleDS(2,(angleDS(2,:) <0)) = angleDS(2,(angleDS(2,:) <0)) + 2*pi;



[C, I] = max(magDS);
hist(I)


hist(magDS(1,:)./magDS(2,:), 50)
hist(magDS(1,:), 50)
hist(magDS(2,:), 50)
hist(magDS{2,1}./magDS{1,1}, 100)
hist(magDS(1,:)./magDS(2,:), 100)

scatter(1:length(angleDS), angleDS(1,:))
scatter(1:length(angleDS), angleDS(2,:))

plot(magDS(1,:), magDS(2,:), '+')
plot(angleDS(1,:), angleDS(2,:), '+')

scatter(1:length(angleDS), angleDS(1,:)./angleDS(2,:))
hist(angleDS(1,:), 30)
hist(angleDS(2,:), 30)
hist(angleDS(1,:)./angleDS(2,:), 30)
hist(angleDS(2,:)./angleDS(1,:), 30)

plot(angleDS(1,:), angleDS(2,:), '+')

plot(magDS(1,:), magDS(2,:), '+')

ginput(1);

Y = get(get(gca,'Children'),'YData');

X = get(get(gca,'Children'),'XData');


h = impoly;

accepted_pos = wait(h);

hold on;

drawPolygon(accepted_pos(:,1), accepted_pos(:,2));

in = inpolygon(X, Y, accepted_pos(:,1), accepted_pos(:,2));

plot(X(in),Y(in),'r+');

cid = 1:183;

chos2 = chos(in);
%chos2 = chos(1,chos2)










