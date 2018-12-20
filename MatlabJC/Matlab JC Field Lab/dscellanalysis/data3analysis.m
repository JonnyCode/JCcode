for i = 2:7
    subplot(2,3,i-1);
    plot(mag{i-1,4}, mag{i,4},'o');
end



for i = 2:7
    subplot(2,3,i-1);
    plot(mag{2,4}, mag{i,3},'o');
end

for i = 4:8
    figure(i)
    subplot(2,3,1);
    plot(mag{4,i}, mag{5,i},'o');
    subplot(2,3,2);
    plot(mag{4,i}, mag{6,i},'o');
    subplot(2,3,3);
    plot(mag{4,i}, mag{7,i},'o');
    subplot(2,3,4);
    plot(mag{5,i}, mag{6,i},'o');
    subplot(2,3,5);
    plot(mag{5,i}, mag{7,i},'o');
    subplot(2,3,6);
    plot(mag{6,i}, mag{7,i},'o');
end


for i = 4:7
    figure(i)
    subplot(2,3,1);
    plot(mag{4,i}, mag{5,i},'o');
    subplot(2,3,2);
    plot(mag{4,i}, mag{6,i},'o');
    subplot(2,3,3);
    plot(mag{4,i}, mag{7,i},'o');
    subplot(2,3,4);
    plot(mag{5,i}, mag{6,i},'o');
    subplot(2,3,5);
    plot(mag{5,i}, mag{7,i},'o');
    subplot(2,3,6);
    plot(mag{6,i}, mag{7,i},'o');
end


for i = 4:8
    figure(i)
    subplot(2,2,1);
    hist(mag{4,i}, 50)
    subplot(2,2,2);
    hist(mag{5,i}, 50)
    subplot(2,2,3);
    hist(mag{6,i}, 50)
    subplot(2,2,4);
    hist(mag{7,i}, 50)
end




plot(mag{4,6}, mag{4,7},'o');
plot(mag{6,4}, mag{7,4},'o');
plot(mag{5,4}, mag{7,4},'o');
plot(mag{6,4}, mag{7,4},'o');
plot(mag{6,4}, mag{7,4},'o');
plot(mag{4,4}, mag{6,4},'o');
plot(mag{4,4}, mag{7,4},'o');
figure(2)

subplot(2,3,1);
plot(mag{4,4}, mag{7,4},'o');
subplot(2,3,2);
plot(mag{4,4}, mag{6,4},'o');    
subplot(2,3,3);
plot(mag{5,4}, mag{7,4},'o');    
subplot(2,3,4);
plot(mag{6,4}, mag{7,4},'o');   
subplot(2,3,5);
plot(mag{6,6}, mag{7,6},'o');    
subplot(2,3,6);
plot(mag{5,6}, mag{7,6},'o');









subplot(2,3,1);
plot(mag{6,6}, mag{7,6},'o');

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

plot(mag{6,6}(1,chos), mag{7,6}(1,chos),'o');



subplot(2,3,2);
plot(mag{4,4}, mag{6,4},'o');    
subplot(2,3,3);
plot(mag{5,4}, mag{7,4},'o');    
subplot(2,3,4);
plot(mag{6,4}, mag{7,4},'o');   
subplot(2,3,5);
plot(mag{6,6}, mag{7,6},'o');    
subplot(2,3,6);
plot(mag{5,6}, mag{7,6},'o');


for i = 1:7
    subplot(7,1,i);
    hist(mag{i,5},50);
end

for i = 1:7
    subplot(7,1,i);
    hist(dsindex{i,4},50);
end

for i = 1:8
    subplot(8,1,i);
    hist(magave(i,:),50);
end

for i = 1:8
    subplot(8,1,i);
    hist(magmax(i,:),50);
end


for i = 2:7
    subplot(2,3,i-1);
    plot(mag{i-1,4}, mag{i,4},'o');
end



for i = 2:7
    subplot(2,3,i-1);
    plot(mag{2,4}, mag{i,3},'o');
end