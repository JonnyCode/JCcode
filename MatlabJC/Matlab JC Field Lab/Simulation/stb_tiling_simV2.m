% starburst tiling simulation
% edit of stb_tiling_sim to include multiple ds cells 

% JC 11/4/2015 
%% parameters
% units are relative to ds dendrite length (1)

% stb center tiling parameters
c_space = .3 ; % distance between stb ceneter nearest neighbors
c_jitter = 0 ; % jitter off of deterministic grid
c_rand = false ; % make random std centers

% stb dendrite parameters
d_number = 20 ; % number of dendrites
d_angle_jitter = 10 ; % degrees of jitter off of evenly spaced dendrites
d_length = .75 ; % length of dendrites
d_sec_length = .01 ; % length of each dendrite integrating section (not important at this point)
d_noise = 0 ; % (degrees) jitter of response (as if d was in different angle than its weighting) 

% stim parameters
space_size = 4 ; % size of stim space
bar_angles = [30:30:360] ; % (degrees) angle of moving bars

% ds rgc parameters
ds_c_x = [2, 1.5, 2, 1.5, 2, 1.5, 2, 1.5] ; % center of ds soma
ds_c_y = [2, 1.5, 2, 1.5, 2, 1.5, 2, 1.5] ; 
ds_nd = [0, 0, 90, 90, 180, 180, 270, 270] ; % (degrees) null direction

%% multiple rounds
for rnd=1;
%% model

% make stb centers (hexagonal grid)
hexfac = sqrt(3)/2 ; % hexagon shift factor
[c_x, c_y] = meshgrid(0:c_space:space_size/hexfac) ; % square matrix
c_x = hexfac * c_x ; % shift proper amount
n = size(c_x,1) ; % size of array
sm = repmat([0 0.5*c_space],[n,ceil(n/2)]) ; % shift matrix
c_y = c_y + sm(:,1:size(c_y,1)) ;

% add stb center jitter
for stb = 1:length(c_x(:)) ; % for each stb
    c_x(stb) = c_x(stb) + normrnd(0,c_jitter) ;
    c_y(stb) = c_y(stb) + normrnd(0,c_jitter) ;
end

% random stb center locations
if c_rand  ;
    c_x = space_size* rand(1,length(c_x(:))) ;
    c_y = space_size* rand(1,length(c_x(:))) ; 
end

% calculate stb dendrite locations
ang = [0:360/d_number:360] ; % angle of dendrite
r=1 ;
for stb = 1:length(c_x(:)) ; % for each stb
    ang2 = ang + rand(1)*360/d_number ; % add global jitter to dendrites of each stb
    for d = 1:d_number ; % for each dendrite
        d_ang(r) = ang2(d) + normrnd(0,d_angle_jitter) ; % angle + jitter of individual dendrite
        d_x(r,:) = c_x(stb) + cosd(d_ang(r))*[0:d_sec_length:d_length] ; 
        d_y(r,:) = c_y(stb) + sind(d_ang(r))*[0:d_sec_length:d_length] ;
        r=r+1 ;
    end
end
  
% calculate stb dendrite weights onto ds
for rgc = 1:length(ds_nd) ;
    for d=1:r-1 ; % for each dendrite 
        syn_fact1(rgc,d) = 1>sqrt((ds_c_x(rgc) - d_x(d,end))^2 + (ds_c_y(rgc) - d_y(d,end))^2) ; % (0,1) end of stb dendrite within range of ds cell?

        SomaA(d) = atan2d(ds_c_y(rgc) - d_y(d,1), ds_c_x(rgc) - d_x(d,1)) ; % angle of soma (-180:+180)
        SomaDa(d) = abs(atan2d(sind(SomaA(d) - ds_nd(rgc)), cosd(SomaA(d) - ds_nd(rgc)))) ; % acute angle between soma and null direction    
        syn_fact2(d) = SomaDa(d)/180 ; % (0-1) gets smaller with smaller acute angle

        dDa(d) = abs(atan2d(sind(d_ang(d) - ds_nd(rgc)), cosd(d_ang(d) - ds_nd(rgc)))) ; % acute angle between dendrite and null direction  
        syn_fact3(d) = dDa(d)/180 ; % (0-1) gets smaller with smaller acute angle

        d_w(rgc,d) = syn_fact1(d) * syn_fact2(d) * syn_fact3(d) ; %(0-1) product of independent weights
    end
end

% calculate dendrite respones
for bar=1:length(bar_angles) ; % for each bar
    for d=1:r-1 ; % for each dendrite
        d_ang_noisy = d_ang(d) + d_noise*rand(1) ; % add noise to dendrite optimal angle
        rA = abs(atan2d(sind(d_ang_noisy - bar_angles(bar)), cosd(d_ang_noisy - bar_angles(bar)))) ; % acute angle between dendrite and bar direction  
        d_r(bar,d) = (rA)/180 ; %(0-1) amplitude of this dendrite to this bar
    end
end

% calculate ds cell response
for rgc = 1:length(ds_nd) ;
    ds_r(rgc,:) = (d_r * d_w(rgc,:)')/sum(syn_fact1(rgc,:)) ; % (0-1) weighted average of dendrites within the rf 
end
        
%% measuring model output
for rgc = 1:length(ds_nd) ;
    ds_vector = PolarVectorAddition([bar_angles', ds_r(rgc,:)]) ; % vector sum
    ds_vectorMag(rnd) = ds_vector(2) ;
    ds_vectorAng(rnd) = ds_vector(1) ;

    [m,i] = max(ds_r(rgc,:)) ;
    ds_peak(rnd) = bar_angles(i) ; % peak of tuning curve

    p1 = interp1(ds_r(1:i),bar_angles(1:i),m/2) ; 
    p2 = interp1(ds_r(i:end),bar_angles(i:end),m/2) ; 
    ds_halfWidth(rnd) = p2-p1 ; % half width of tuning curve in degrees
end


% figures

%stb_plot = [1000,10000,30000,30300,40000] ;
stb_plot = [1000] ;

figure

subplot(2,2,1)
viscircles([ds_c_x,ds_c_y],1) ;
hold on
plot(c_x(:),c_y(:),'.k')
for a= 1:length(stb_plot) ;
    for d = [1:d_number]+stb_plot(a) ;
        plot(d_x(d,:),d_y(d,:),'color',[syn_fact1(d), 0, 0])
    end
end
title('Syn fact 1')

subplot(2,2,2)
viscircles([ds_c_x,ds_c_y],1) ;
hold on
plot(c_x(:),c_y(:),'.k')
for a= 1:length(stb_plot) ;
    for d = [1:d_number]+stb_plot(a) ;
        plot(d_x(d,:),d_y(d,:),'color',[syn_fact2(d), 0, 0])
    end
end
title('Syn fact 2')
    
subplot(2,2,3)
viscircles([ds_c_x,ds_c_y],1) ;
hold on
plot(c_x(:),c_y(:),'.k')
for a= 1:length(stb_plot) ;
    for d = [1:d_number]+stb_plot(a) ;
        plot(d_x(d,:),d_y(d,:),'color',[syn_fact3(d), 0, 0])
    end
end
title('Syn fact 3')

subplot(2,2,4)
viscircles([ds_c_x,ds_c_y],1) ;
hold on
plot(c_x(:),c_y(:),'.k')
for a= 1:length(stb_plot) ;
    for d = [1:d_number]+stb_plot(a) ;
        plot(d_x(d,:),d_y(d,:),'color',[d_w(d), 0, 0])
    end
end
title('dendrite weights')

figure
subplot(1,3,1)
plot(bar_angles, ds_r)
hold on
xlabel('bar angle')
ylabel('starburst input')

subplot(1,3,2)
hist(d_w,[0:.01:1])
hold on
xlabel('dendrite weights')
ylabel('number observations')

subplot(1,3,3)
hist(d_ang,[0:360])
hold on
xlabel('dendrite angles')
ylabel('number observations')

figure(4)
subplot(1,2,1)
plot(bar_angles, ds_r)
hold on

subplot(1,2,2)
polar([0,ds_vectorAng(rnd)*pi/180],[0,ds_vectorMag(rnd)])
hold on

end

%% histograms

[ds_vectorAng_histY, ds_vectorAng_histX] = hist(ds_vectorAng,[0:5:360]) ;
[ds_vectorMag_histY, ds_vectorMag_histX] = hist(ds_vectorMag,[0:.01:1]) ;
[ds_peak_histY, ds_peak_histX] = hist(ds_peak,[0:5:360]) ;
[ds_halfWidth_histY, ds_halfWidth_histX] = hist(ds_halfWidth,[0:5:360]) ;

% figures
figure
subplot(1,4,1)
plot(ds_vectorAng_histX, ds_vectorAng_histY)

subplot(1,4,2)
plot(ds_vectorMag_histX, ds_vectorMag_histY)

subplot(1,4,3)
plot(ds_peak_histX, ds_peak_histY)

subplot(1,4,4)
plot(ds_halfWidth_histX, ds_halfWidth_histY)


%% for future code maybe...   

% % make moving bar stimulus
% for bar=1:length(bar_angles) ; % for each bar
%     for t = 1:length(time) ; % for each time point
%         for a = 1:space_size/pix_size ; % for each pixel X
%             for b = 1:space_size/pix_size ; % for each pixel Y
%                 temp(a,b) = a*cosd(bar_angles(bar)) + b*sind(bar_angles(bar)) - time(t)/bar_speed  ;
%             end
%         end
%         stim{bar}(:,:,t) = (temp>=-bar_width/2 & temp<=bar_width/2) ;
%     end
% end

% % calculate response of dendrite to moving bar
% for d=1:r-1 ; % for each dendrite
%     for bar=1:length(bar_angles) ; % for each bar
%         for t = 1:length(time) ; % for each time vector in response model
%             d_stim(t,:) = interp2([0:pix_size:space_size],[0:pix_size:space_size], stim{bar}(:,:,t), d_x(d,:), d_y(d,:)) ;
%         end % FIXXXXXX
%     end
% end


    