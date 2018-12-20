%       This script will help you better understand the cross correlation
% function.  You will also need the Decimate.mat, and RGC_example_data files. 

% Jon Cafaro 3/1/07

% What is explained in this tutorial?
%   The purpose of this tutorial is to explain the cross correlation function, 
% and give an example of its use.  When you are done hopefully you will have 
% a clearer understanding of what aspects of two signals may be reflected 
% in a cross correlation and how that may be calculated.

% What is the Cross Correlation function?
%       The Cross Correlation function is a measure of the correlation
% between two signals as they are time shifted against each other.  This 
% measure can help us to determine how similar two signals are from each 
% other and in what ways they may differ.  

%       Lets start with the correlation of two identical signals at a given time 
% shift of zero.  This is actually called an autocorrelation because it is the same
% signal but this will start to give us an idea for the types of numbers we
% may see when two signals are similar.


xx = [4,3,2,1,2,3,3,3,2,1,4,3,2,1,2,3,3,3,2,1];      % this is a really simple descreet signal (you could actually do these calculations on paper)
tt = [1:length(xx)];                                 % we can say each point was taken at each time point tt
figure, plot(tt,xx)                                 % take a look

yy = xx ;                                             % We'll make the second signal the same as the first 

% Correlation = (<xy>-<x><y>)/sqrt(var(x)*var(y))
% this is represented below

Covariance = mean(xx.*yy)-mean(xx)*mean(yy)  % this takes the covariance of the two signals


normalization = sqrt((std(xx,1))^2*(std(yy,1))^2) % this factor will normalize our correlation 
                                                % NOTE: the standard deviation calculated here is slightly 
                                                % different than what matlab will do by default but this is 
                                                % a necessity because we are dealing with such a small signal
                                              
                                                
correlation = Covariance/normalization         % this is the correlation given no time shift (lag = 0)

% NOTICE - the correlation is equal to 1.  The highest cross correlation
% coefficient we can have is 1 (When two signals are the same).
% NOTE- this is not the case if we are not using a normalization factor 
% (ie. if we're not looking at the coeficient)

% Now lets start to shift these vectors in time against each other and see
% what happens to the correlation.
                 
yy_shift = [zeros(1,length(xx)),yy,zeros(1,length(xx))]             % One way to slide the signals against each other is 
                                                                    % to add zeros so they cannot slide off each other.    
xx_shift = zeros(1,length(yy_shift))

for time_shift = 1:length(xx)+length(yy)                                     % we will shift the signal xx all the way over signal yy 
yy_shift = [zeros(1,length(xx)),yy,zeros(1,length(xx))]
xx_shift = zeros(1,length(yy_shift))

xx_shift(1,time_shift:time_shift+length(xx)-1)=xx                    


Covariance = mean(xx_shift.*yy_shift)-mean(xx_shift)*mean(yy_shift)          % this is as above
normalization = sqrt((std(xx_shift,1))^2*(std(yy_shift,1))^2)
correlation(time_shift) = Covariance/normalization
end

tt2=1:length(correlation)                                % this will create a time vector that represents where xx and yy are in respect to each other
tt2=tt2-(length(tt2)+1)/2             
figure, subplot(2,1,1)                                   % plot cross corr and signal
plot(tt2,correlation)
title('cross Correlation')
xlabel('lag(ms)')
ylabel('cross correlation coefficient')
subplot(2,1,2)
plot(tt,xx)
title('signal')
xlabel('time')
ylabel('amplitude')

% Look at the cross correlation function and the signal that it reflects.
% NOTICE:   1) It is symetric (although only because it is an autocorrelation- a cross corr may not be symetric )
%           2) At lag 0 it has a peak of 1, as we noted earlier
%           3) There are other small peaks - can you tell why these peaks exist?
%           4) There are negative correlations - so the cross correlation can also
%              reflect events in the signals that are anticorrelated (ie. when one is up
%              the other is down)
%           5) The cross correlation has a general triangular shape - this
%              is a result of the zeros we added and we will see that some of
%              the matlab functions create similar effects.

% How does Matlab calculate cross correlation?
%    Matlab uses a function called "xcorr", but unlike our example above this
% does not simply calculate the correlation at various time shifts.
% Instead the heart of the xcorr function is a fourier transform.

% What does the fourier transform have to do with a cross correlation?
%   The correlation of two continuous functions is:
%           The Integral of f(t)h(t+T)dt
%   This is very similar to a convolution, which is:
%           The Integral of f(t)h(t-T)dt
% You'll notice the only difference beween the convolution and correlation
% is that the convolution flips the second signal before it shifts it (-T)
% while the correlation does not (+T).  To solve a convolution we can
% take the fourier transorms of each signal, multiply them and then take
% the inverse of this value, as such:
%       Convolution = Inv(fft(f(t))fft(h(t-T)) = Inv(f(w)h(w))
% We can do a very similar calculation to solve the correlation but because
% the correlation does not flip the second signal we must use the complex
% conjugate of h(w).  Matlab uses the following peice of code to calculate
% the cross correlation:
%       Cross Correlation = Inv(fft(f(t))*conj(fft(h(t+T)))) = Inv(f(w)conj(h(w)))
        
% Lets use a more realistic example of signals on which we may want to run
% a cross correlation.  

load RGC_example_data                               % loads examples of spike spike trains from retinal ganglion cells 

figure,
subplot(2,1,1)
plot(ON_CELL_545tutorial)
axis([0 85000 0 1.5])

hold on,plot(Offt_CELL_545tutorial,'r') 
hold on,plot(Offs_CELL_545tutorial,'g') 
title('spikes response to light fluctuations')
xlabel('time (0.1ms)')
ylabel('spike')

subplot(2,1,2)
plot(ON_CELL_545tutorial)
axis([60000 65000 0 1.5])

hold on,plot(Offt_CELL_545tutorial,'r') 
hold on,plot(Offs_CELL_545tutorial,'g') 
title('spikes response to light fluctuations')
xlabel('time (0.1ms)')
ylabel('spike')
legend('On Cell', 'Offt Cell', 'Offs Cell')

%   The preceding figure shows the spike times of 3 different kinds of
% retinal ganglion cells in response to a fluctating light stimuli.  The On
% cells will respond to the onset of light a step while the off cells will
% respond to the offset.  Furthermore, the offs (off-sustained) respond with 
% longer lasting spike trains than the offt (off transient) and seem to 
% respond consistently earlier.    

%   Let's see what there cross correlations will look like.


CrsCorrs2(1,:) = xcorr(Offt_CELL_545tutorial,ON_CELL_545tutorial,'coef');     % this calls the matlab function xcorr and uses the 'coeff' option
CrsCorrs2(2,:) = xcorr(Offs_CELL_545tutorial,ON_CELL_545tutorial,'coef');
CrsCorrs2(3,:) = xcorr(Offs_CELL_545tutorial,Offt_CELL_545tutorial,'coef');
CrsCorrs2(4,:) = xcorr(Offt_CELL_545tutorial,Offs_CELL_545tutorial,'coef');

for a=1:4
CrsCorrs(a,:) = DecimateWave(CrsCorrs2(a,:), 10);                         % this averages some of the cross corelation values to produce a smoother graph
end

tt = 1:length(CrsCorrs);                                                % prepares a time axis
tt=tt-length(CrsCorrs)/2;
figure, subplot(4,1,1)                                                  % plots
plot(tt,CrsCorrs(1,:));
subplot(4,1,2) 
plot(tt,CrsCorrs(2,:),'r');
subplot(4,1,3) 
plot(tt,CrsCorrs(3,:),'g');
subplot(4,1,4) 
plot(tt,CrsCorrs(4,:),'k');


%    Look at the four cross correlations.  NOTICE: 1) That they all have that
% triagular shape, this is a result of zero padding that the xcorr function
% performs, similar to the zero padding that we did in our simple example.
% This zero padding causes the the correlations to become smaller and
% smaller with greater shifts (further from lag=0) because the signals are
% being correlated with more and more zeros.  2) Also notice that unlike
% our previous example there are no negative correlations here.  This is
% because the xcorr function does not subtract the mean as we did in our
% previous example (remmeber we used <xy>-<x><y>) and so will not produce 
% negative correlation values unless some of the signal has negative values.
% Both of these problems can be fixed if we use 'xcov' instead of 'xcorr' 
% but these are details and xcorr has still preserved the important 
% parts of our cross correlations. 


% The peaks of the cross correlation functions can be important to focus
% on.  So let zoom in on the peaks.

subplot(4,1,1)
axis([-400 400 0 .014])
subplot(4,1,2)
axis([-400 400 0 .014])
subplot(4,1,3)
axis([-400 400 0 .014])
subplot(4,1,4)
axis([-400 400 0 .014])

% Notice - 1) The top two graphs (blue and red) have dips near the zero lag.  This is because
% we are comparing an on cell with an off cell in both of these cross
% correlations.  Therefor, we might expect to see that near the time an off
% cell is spiking (lag=0) our on cell would have a good chance of not
% spiking, producing a low correlation value.  2) The third graph has a
% peak rather than a dip near the zero lag.  This is also expected because 
% they are both off cells.  You'll notice that this peak is slightly off zero, 
% in fact it seems to occur about 50ms before.  This indicates that the offs cell was
% likely to spike about 50 ms before the offt cell.  3)  Compare the third
% and fourth graphs.  Notice they are both comparing an offt to an offs
% cell but they are comparing them in different orders.  These graphs
% should have the same peak but opposite lags.  In fact, these graphs are
% mirror images of each other only reflected about the zero lag axis.  Try using 
% the 'flip' command on one of the traces and you can varify this fact for yourself. 
% Disclaimer- they will not be exact mirror images because of effects of
% the Decimate command, but they should be very close.  4)  If you compare
% the peak correlation values of these comparisons you will see that offt
% vs offs have a higher peak than on vs off cells.  This would indicated
% that these cells are more similar to eachother in their spike timing than
% the other pairs.

% We can also get information from the width of the peak and antipeaks of
% the cross correlation functions, but this is a more convoluted process
% and requires use of the auttocorrelation as well.  So it will not be
% discussed.

% Now that we have some sense for what information a cross correlation 
% function can provide lets go back to see what the cross correlations will 
% look like if we do not use xcorr but instead simple plot what I
% previously refered to as the "heart of the xcorr function".

% As mentioned the xcorr function can be described by:

a=fft(ON_CELL_545tutorial);              % the fourier transform of signal a
b=fft(Offt_CELL_545tutorial);            % the fourier transform of signal b


raw_xcor = ifft(a.*conj(b));             % the inverse of signal a times the complex conjugate of signal b
                                        % NOTE the little dot after the a,what does this do? try testing it with simple matices if you don't know 
                                       

raw_xcor= DecimateWave(raw_xcor, 10);   % again just decimating to smooth things out, this wouldn't be such a problem if we had more example data to average

% check put the result
time=1:length(raw_xcor);
figure, plot(time,raw_xcor);

% Things are different! Fistly, that is because the fourier transform has not
% centered our peak, in fact its at the end of our plot.  Let fix this

raw_xcor2 = [raw_xcor((length(raw_xcor)/2+1):end),raw_xcor(1:(length(raw_xcor)/2))];

time=time-(length(time)+1)/2 ;
figure, plot(time,raw_xcor2)

% Okay now this looks better, but it still looks different then what xcorr
% produced.  NOTICE-1) it does not have that triagular shape that xcorr
% produced. This is because there is no zero padding in this case.  2) It
% shows cross correlation values above 1.  This is because we have not 
% normalized (ie. these are not coeficients).  

% lets compare the peaks of the cross corr we get with xcorr and the one we
% get with our fourier transforms.

figure, subplot(1,2,1)
plot(time,raw_xcor2)
subplot(1,2,2)
plot(tt,CrsCorrs(1,:));

subplot(1,2,1)
axis([-400 400 0 3.5])
title('fourier raw version')
xlabel('lag(ms)')
ylabel('correlation value')
subplot(1,2,2)
axis([-400 400 0 .01])
title('matlabs xcorr')
xlabel('lag(ms)')
ylabel('correlation value')

% not bad, the major features seem in tact.  You may find, as I have, that
% using different methods of calculating the cross corelation will produce
% slightly different plots but that the important parts of those plots are
% maintained.  So hopefully you now have a better understanding of what 
% information can be learned from the cross correlation function and how 
% the function can be calculated. 
