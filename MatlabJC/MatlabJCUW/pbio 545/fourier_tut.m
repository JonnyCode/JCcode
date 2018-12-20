% FourierTutorial.m
%****************************************************************************
% Fourier analysis tutorial for NeuBeh/PBIO 545, Winter 2003
%
% Created 1/03 Mike Shadlen
% Revisions: 1/05 update
% Part I revised by Adrienne 1/06
%
%****************************************************************************

% The Fourier transform is probably the most important transform in applied
% math. It takes a function, typically of time or space, and expresses it
% as a function of frequency. The goals of this tutorial are to make you
% comfortable with what a Fourier transform is, how to compute one, and why 
% it is useful. We will focus on two topics: (1) definition and properties of 
% the Fourier transform, (2) convolution, and (3) the connection between the two. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part I. The discrete Fourier transform

% A transform is (usually) a reversible mapping from one set of values to
% another set of values. The Fourier transform takes a list of intensities 
% at points in time (or space) and transform them to a list of "intensities" 
% at different frequencies.
% Thinking about a signal as a collection of time intensities, each
% associated with its own time point, is simple enough. Thinking about a
% signal as a list of intensities, each associated with its own sinusoidal
% frequency, is also convenient. Perhaps because of the way our auditory system 
% is organized, we find thinking in frequency quite natural and intuitive.
% You do it every day when you adjust the
% bass and treble on your stereo -- or better yet, if you adjust knobs on
% an equalizer. You can imagine that describing a signal in terms of its
% frequencies can be convenient in certain settings. That's what the
% Fourier transform is for. 
% 
% We should emphasize from the start that taking a Fourier transform (FT, for
% short) does not change the signal; it just represents it differently. To
% be concrete, suppose there is some electrical signal that we want to
% describe over the phone to a friend. We could list the voltages at each
% time or we could list the amplitudes and phases of each of the
% frequencies that comprise the signal. It's the same signal, just a
% different list of numbers. Sometimes it is much faster and simpler to use
% one type of description than the other.
% 
% Often, when we describe signals, we don't think of them as lists of
% numbers but rather as functions of time. Of course, the function can be
% viewed as a way to generate a list of values at each time point. We
% can also use a function to describe the list of "intensities" at each
% frequency in a Fourier transform. The two functions are called Fourier
% transform pairs. For example, suppose s(t)  is some function of time (it
% has a value at each time point). And suppose that you could describe
% the same signal by "intensities" of frequency (w) using a function, S(w).
% Then, we say that s(t) and S(w) are a Fourier transform pair. Or we say
% that S(w) is the Fourier transform of s(t), or that s(t) is the inverse
% Fourier transform of S(w).  When we started by saying that a transform is
% a reversible mapping, what this means is that S(w) is the unique Fourier
% transform of s(t) and s(t) is the unique inverse transform of S(w).

% Our goal here is both to help you to get comfortable with the meaning of
% the Fourier transform, and to get you used to using the numerical
% function FFT (fast Fourier transform) in Matlab. If you get everything
% about what that matlab function returns to you, you'll pretty much have 
% mastered what you need to know.

% Before we look at FTs, let's get comfortable with a few illustrative
% functions. Since we'll be doing everything numerically in this tutorial,
% we will focus here on discrete functions. What is meant by this is
% that we are going to make the argument of the function, say time, take 
% discrete values. We're used to thinking of time t occupying some continuous 
% interval on the real number line, but here we'll be thinking of time as 
% integers, multiplied by some constant or scaling factor. In fact this is 
% overwhelmingly going to be the case whenever you are thinking about real 
% data: the values are always sampled at some particular finite sampling rate. 

clear all

% Define a discrete time axis. The spacing is dt
dt = .01;      % our units are seconds. So each integer time point is scaled to represent dt sec
tmax = 1
t = [0:dt:tmax-dt]';
% Notice that we are sampling time at a rate of 1/dt. This will be
% important later.
samplingRate = 1/dt    % in Hz

% Think of a discrete function as a series of weights at each of the points
% in t.  Let's look at some functions this way. Here's a gaussian profile
% centered at .5
f1 = exp(-.5*(t-.5).^2 / (.05)^2);
figure(1),clf,hold on
h = plot(t,f1)
xlabel('Time [s]')
ylabel('Amplitude')
% Pause here. The function is drawn continously: there are no gaps. That's
% because matlab drew a line, or interpolated, between the values at the discrete time
% points. In reality, what we know about the function is just a bunch of heights at each of
% the time points. We can use the stem function to highlight this fact.

stem(t,f1,'filled');

% and we can toggle the smooth curve off
set(h,'Visible','off')

% Now... that's a discrete function. Each of the points can be thought of as
% a weight given to the basis set of discrete time points. We represent the
% function as a list of these weights. Indeed, f1 is a vector.
size(f1)

% Let's remember this function by calling it something mnemonic
gaus1 = f1;

% Let's look at a few more functions
% pulse
f2 = zeros(size(t));
f2(t>.45 & t<= .55) = 1;
figure(1),hold off
stem(t,f2,'filled')
xlabel('Time [s]')
ylabel('Amplitude')
pulse1 = f2;

% sinewaves
f3 = sin(2*pi*1*(t-.5));
sin1 = f3;
stem(t,f3,'filled')
%
f4 = cos(2*pi*10*(t-.5));
cos10 = f4;
stem(t,f4,'filled')
% This looks a lot better interpolated!
plot(t,f4)

% the sum (I'm not giving this a function name)
sos = .5*f3 + f4;
plot(t,sos)

% exponential
tau = .1
f5 = exp(-(t-.5)/tau);
f5(t<.5) = 0;
exp1 = f5;
stem(t,f5,'filled')

% discrete delta function -- take special note of this one. 
f6 = zeros(size(t));
f6(min(find(t>=.5))) = 1;
stem(t,f6,'filled')
delt1 = f6;

% comb(.1t) function. This is a sampling function or strobe function that
% turns on at a rate of 10 times per second (10 Hz). But it is also a sum
% of discrete delta functions.

f7 = mod(t,.1)==0;
comb10 = f7;
stem(t,comb10,'filled')

% Reflection point. You have yet to see a Fourier transform. You are
% supposed to simply marvel at the idea of discreteness. And here is the
% one thing extra that you need to think about. All of the functions above
% can be thought of as sums of discrete delta functions with appropriately chosen
% weights. The most obvious example is the comb function, as we have
% already noted. But all the functions are just values at discrete time
% points. So all can be thought of as sums of delta functions. In fact, a discrete
% function is exactly this: the multiplication of a continuous function by
% comb(d*t), where d is the gap between samples (i.e., the reciprocal of the
% sampling rate). We will return to this point later.
%
% Another important point to make before moving on. Let's recall the idea of a
% BASIS SET from linear algebra.  As you know, a coordinate frame is a set of 
% vectors which span a given space: e.g. in the plane, we have an x axis, and 
% a y axis. So to say that a point in the plane is (3,2) means that our vector 
% contains 3 units along the x axis, and 2 units along the y axis. 
% A functional basis set is just a coordinate set for the space of
% functions. What does this mean?  Let's take some arbitrary data sequence, a bunch of
% values from time zero to time T. Since it is discretized, we have N points 
% in our data: a value for at each time, 0, dt, 2dt, 3dt,.. (N-1)dt.
% Let's now think about this data as an N-dimensional vector, and plot it
% as a single point in an N-dimensional space.  What is the coordinate
% system in which we have plotted our function? The unit vectors are the set 
% of delta functions at each time, t = i dt, for i ranging from 0 to N-1. So axis one
% means how strong was the function at time t = 0, axis 2 is how strong
% the function was at time t = dt, and so on.  If we write these as vectors
% they will look just like Euclidean coordinate axes: [1 0 0 0 0 0 ..], [0 1 0 0
% 0 ...], etc. 
% Remember that an important property of many useful coordinate systems is
% that the axes are orthogonal: the x vector has exactly zero projection onto 
% the y vector; the dot product of the x vector with the y vector is zero.
% Similarly for our delta function coordinate system. The
% delta function at time t is independent of the delta function at time t'
% if t is not the same as t'; the dot product of these two functions is zero.
% As we said already, any function can be written in this "time" basis by summing 
% up weighted delta functions.

% Now, let's move onto Fourier transforms. We had better start by defining 
% it.  As we've said, calling it a transform already gives away something 
% important. The function or data we're describing is not changed; it's just expressed 
% differently. Just like when you rotate something --you don't change the structure 
% of the thing you rotated. 
%
% So now, instead of thinking of the function as a list of weights at the
% discrete time samples, we would now like to think of it as a list of "weights" on
% sinusoids. In other words, since the time representation was just one 
% choice of coordinate systems, we can just as well represent our data function, 
% or N-dimensional point, in some other coordinate system. The Fourier transform 
% is a method for representing our data in another set of axes, where the axes are  
% sines and cosines rather than delta functions. 

% This is done in a way that will need a little unpacking. Let's look at the definition
% of the Fourier transform F(w) of a continuous function f(t):
%
% F(w) = (1/sqrt(2 pi)) integral dt f(t) exp( - i w t). 
%
% (check out the wikipedia entry to see this in proper math type!)

% The integral here is the equivalent for a continuous function of a dot product or 
% projection. So the idea is that we project our data function against the basis 
% function exp(-i w t), where w is frequency; one basis function for every value of w. 
% So how is this sines and cosines? Remember back to the olden days when you first met 
% complex numbers. Remember Euler's identity:
%
% exp( i x) = cos(x) + i sin(x). 
%
% Let's have a look at some examples:

exp(i*0.5)

exp(i*pi)

exp(i*-pi)

exp(i*pi/2)

% 
% All of these are complex numbers, but with different real and imaginary
% parts. (If the above didn't give you back a complex number, try 

clear i

% in case you have redefined "i" somehow as a real number.)

% I included some particular cases where the real or the imaginary
% part is zero. So what is this i thing? It's telling you
% that this is really a two-component number, with two "axes", the "real"
% axis and the "imaginary" axis, where the length of the number exp(i x) along the 
% real axis is cos(x), and its length along the imaginary axis is sin(x).
% We could write it in vector notation as (cos(x), sin(x)) and plot it as a
% point on a plane. The length of the vector is 1, since 
%
% cos^2(x) + sin^2(x) = 1. 
%
% If we want a longer vector, we scale the whole thing by some amplitude
% A: A exp( i x). The value of x is the phase in radians.
%
% So: every complex number is really 2-dimensional, and every  exp(- i w t) 
% for a given w will give us two components, one real, corresponding to the
% cosine at frequency w, and one imaginary, corresponding to the sine.
% This is just an elegant way to capture the issue that for every frequency, we need two
% values. Why? Because every frequency component is specified both by an
% AMPLITUDE and a PHASE. Using both a sine and a cosine component allows us
% to construct any arbitrary phase. Remember another fact from your dim
% dark past: the sin identities. Do you remember:
%
% sin (x+y) = sin(x) cos(y) + sin(y) cos(x)?
%
% so let's say x is t and y is some phase p. 
% then a sin (t + p) = a * cos(p)*sin(t) + a * sin(p)*cos(t), where a*cos(p) 
% and a*sin(p) now become the coefficients of the sine and cosine functions. Again, 
% the total amplitude is a * sqrt(sin(p)**2 + cos(p)**2) = a!  So representing phase 
% as well as overall amplitude is why you need both a sine and a cosine at every 
% frequency. Using complex numbers for each frequency is just a way to do
% this neatly and elegantly.

% Let's see that in action. Let's look at the real and imaginary parts of the function 
% exp(i w t), for w = 2 pi.  (this frequency is angular frequency! here we
% will generally take out the 2 pi factor, w = 2 pi f, and work in hertz..)

figure(2), clf
subplot(2,1,1); plot(t,real(exp(i*2*pi*t))); xlabel('t'); title('Real part of exp(iwt)')
subplot(2,1,2); plot(t,imag(exp(i*2*pi*t))); xlabel('t'); title('Imaginary part of exp(iwt)')

% Let's also look at a fun demo showing how you shift the phase of a sine
% or cosine function by multiplying by exp(i phi), where phi is the phase.
% Below the graph you'll see the polar representation of exp(i phi).
%
figure(4), clf
clear ax
ax(1) = subplot(3,1,1);
ax(2) = subplot(2,3,5);
axis square
phi = [0:.01:2*pi]';  % take phase shifts going right around the circle.
polar(phi,ones(size(phi)));
polar([0 0],[0,1])
set(gcf,'CurrentAxes',ax(1))
plot(t,cos(2*pi*t))
set(gcf,'CurrentAxes',ax(2))
nsteps = 100;         % change this to alter the speed of the demo
for j = 1:nsteps
    phi = j* 2*pi/nsteps;  % we are incrementing the phase shift
    set(gcf,'CurrentAxes',ax(1))
    plot(t,real(exp(i*(2*pi*t + phi))),'r','LineWidth',3)
    set(gcf,'CurrentAxes',ax(2))
    polar([0 phi],[0 1],'r'),hold on
    polar([0 0 phi],[0 cos(phi) 1],'k--'), hold off
    % pause
    drawnow
end

% So what's the upshot: by projecting onto exp(i w t), we are actually 
% projecting against both a sine and a cosine function of frequency w.  
% Remember that a sine, or any other function, is just a linear sum of the 
% delta functions so this is just a linear transformation. 
% Even better, since sines and cosines of different frequencies are orthogonal, 
% the sines and cosines form an orthogonal basis set. The sines and cosines are not the only
% functions that we could choose for our alternative basis set, but they
% are useful so often that the Fourier transform is one of the most
% commonly used of all possibilities. When we study dimensionality
% reduction, we will talk about some other alternatives.

% Let's look at some examples.

% Let's start with a cosine function. What do you expect to see? Let's bung
% the cosine function into matlab's function fft and plot. I will use a
% frequency f of 10 Hz. Worth mentioning now that in the above we have used the 
% symbol w for omega, which is the angular frequency. Frequency in Hz is related to
% angular frequency by a factor of 2 pi, w = 2 pi f, which is often useful to factor
% out. 

y = cos(2*pi*10*t);
figure(2), clf
subplot(1,2,1)
plot(t,y)
xlabel('Time')
subplot(1,2,2)
fcos = fft(y);
plot(fcos); ylabel('Fourier transform of cosine function')

% What the *$&#^%?  (pardon my Australian.)
% Let's see what you have plotted.

whos fcos

% Ah ha. The fft of y is, of course, a complex array, so has a real and imaginary part.
% When you use plot, matlab plots the real part vs the imaginary part.
% So let's plot each of the two parts separately.

subplot(1,2,1); 
plot(real(fcos));
ylabel('Real part of fft(cos t)');
subplot(1,2,2); 
plot(imag(fcos));
ylabel('Imaginary part of ft(cos t)');

% Now we're getting somewhere. The real part has two big spikes; the imaginary 
% part is tiny, really just zero except for numerical noise.
% Now, why and where are those two big spikes?
% We put in a 100-dimensional purely real data vector, and got out a
% 100-dimensional complex vector, meaning 200 values. That seems to be
% twice as many as we need. We'll talk more about this later. 
% We have obtained a function of frequency, but what are the labels on the 
% frequency axis? 

% Let's recall yet another high school fact: how do you
% write a cosine in terms of complex numbers? It's just a rewriting of the
% definition (known as the Euler identity) we used above. 

% Remember: 
%
%   cos(a) = [ exp(ia) + exp(-ia) ]/2 .
% 
% Similarly, 
%
%   sin(a) = [ exp(ia) - exp(-ia) ] / 2i.
%
% Our transform scans over w (or f), looking for nonzero values. Under the
% integral sign, we have
% 
% exp(i w t) [ exp(i 2 pi 10 t) + exp( -i 2 pi 10 t)] = 
%               exp(i (w- 2 pi 10) t)) + exp(i(w + 2 pi 10) t). 
%
% The integral of an imaginary exponential turns out to be zero unless what is 
% in the exponent exactly cancels, giving you exp(0) = 1. This happens for 
% the first term at w = 2 pi 10, and for the second term at w = - 2 pi 10.

% So, we expect two peaks for our transform, at f = 10 and at f = -10. 

% Why do the peaks show up in the output array at index 10 and 90? It turns
% out that the negative frequencies are arranged after the positive ones. 
% The ordering runs 0 through to max frequency W, then from -W to -min
% frequency. From now on, we'll use the matlab function "fftshift" to rearrange 
% them in the usual increasing order. 
%
% Let's have a look at the transform of a sine wave as well:
% (note here I'm switching to a plotting form more appropriate for this
% kind of spiky function, and I'm scaling the y axes of the real and imaginary
% parts by the same amount so we don't see all that 10^-13 stuff)

clf;
subplot(2,2,1); 
stem(fftshift(real(fcos)),'filled');
ylabel('Real part of ft(cos)');
subplot(2,2,2); 
stem(fftshift(imag(fcos)),'filled');
ylabel('Imaginary part of ft(cos)');

y = sin(2*pi*10*t);
subplot(2,2,3)
fsin = fft(y);
subplot(2,2,3); 
stem(fftshift(real(fsin)),'filled');
ylabel('Real part of ft(sin)');
subplot(2,2,4); 
stem(fftshift(imag(fsin)),'filled');
ylabel('Imaginary part of ft(sin)');

% Just as we should now expect: the large values are in the imaginary part,
% because of that i in the definition of sine, and include a positive and 
% a negative peak, because of the sign difference between the exponents.

% Did you expect that the +f component would come out with a negative sign? Why?

% This is a good moment to think what exactly the minimum and maximum 
% frequencies are. 
% The lowest frequency is determined by the length of the sample. If the
% sample is of length T, then the lowest frequency that can be resolved is
% 1/T.  How about the highest frequency?  If the sampling rate is dt, then
% you cannot specify any frequency higher than 1/2dt. That is because you
% need at least 2 points to capture an oscillation, one on the up part and
% one on the down part. This frequency is called the NYQUIST frequency.
% The Nyquist frequency says that you must sample at least twice as often as 
% the shortest wavelength in your data.  In class we'll try to look at what happens 
% if your signal contains frequencies higher than your sampling rate.

% OK, so let's define our frequency axis.

nyq = samplingRate/2;
dw = 1/tmax;
fax = -nyq: dw: nyq-dw;

% and now finally plot, 

clf;
subplot(2,2,1); 
stem(fax,fftshift(real(fcos)),'filled');
set(gca,'YLim',[-50 50])
ylabel('Real part of ft(cos wt)');
xlabel('Frequency')
subplot(2,2,2); 
stem(fax,fftshift(imag(fcos)),'filled');
set(gca,'YLim',[-50 50])
ylabel('Imaginary part of ft(cos wt)');
xlabel('Frequency')

y = sin(2*pi*10*t);
subplot(2,2,3)
fsin = fft(y);
subplot(2,2,3); 
stem(fax,fftshift(real(fsin)),'filled');
set(gca,'YLim',[-50 50])
ylabel('Real part of ft(sin wt)');
xlabel('Frequency')
subplot(2,2,4); 
stem(fax,fftshift(imag(fsin)),'filled');
set(gca,'YLim',[-50 50])
ylabel('Imaginary part of ft(sin wt)');
xlabel('Frequency')

% OK, let's stretch our wings a bit and try some other functions.

% A sum of 3 sinusoids!

y = sin(2*pi*3*t) + .33 * sin(2*pi*9*t) + .2 * sin(2.*pi*15*t); 
figure(3), clf
subplot(2,2,1)
plot(t,y)
xlabel('Time')
subplot(2,2,2)
stem(fax,fftshift(imag(fft(y))),'filled')
xlabel('Frequency')

% Notice that the time function is almost square wave like. Its definition
% makes it crystal clear that it is a sum of three sinusoids. This is
% pretty easy to see in the Fourier transform. We only plotted the
% imaginary part as we now know that sines will only give us nonzero imaginary
% part.

% Now let's look at each of the functions we've talked about already. For
% ease of viewing, I will plot some of the time functions with lines that
% interpolate between the (time or frequency) points. For some functions, I
% use the stem command to emphasize the discreteness.  First, let's just
% plot the overall amplitude (the absolute value of the complex number).

figure(3), clf
k = 1;
ax(k) = subplot(7,2,k); k=k+1;
hg(1) = plot(t,f1);
ax(k) = subplot(7,2,k); k=k+1;
hg(2) = plot(fax,fftshift(abs(fft(f1))));
% stem(fax,fftshift(abs(fft(f1))),'filled')
ax(k) = subplot(7,2,k); k=k+1;
hg(3) = plot(t,f2);
ax(k) = subplot(7,2,k); k=k+1;
hg(4) = plot(fax,fftshift(abs(fft(f2))));
% stem(fax,fftshift(abs(fft(f2))),'filled')
ax(k) = subplot(7,2,k); k=k+1;
hg(5) = plot(t,f3);
ax(k) = subplot(7,2,k); k=k+1;
stem(fax,fftshift(abs(fft(f3))),'filled')
ax(k) = subplot(7,2,k); k=k+1;
hg(6) = plot(t,f4);
ax(k) = subplot(7,2,k); k=k+1;
stem(fax,fftshift(abs(fft(f4))),'filled')
ax(k) = subplot(7,2,k); k=k+1;
hg(7) = plot(t,f5);
ax(k) = subplot(7,2,k); k=k+1;
hg(8) = plot(fax,fftshift(abs(fft(f5))));
ax(k) = subplot(7,2,k); k=k+1;
stem(t,f6,'filled')
ax(k) = subplot(7,2,k); k=k+1;
stem(fax,fftshift(abs(fft(f6))),'filled')
ax(k) = subplot(7,2,k); k=k+1;
stem(t,f7,'filled')
xtick = get(gca,'XTick');
xticklabel = get(gca,'XTickLabel');
xlabel('Time')
ax(k) = subplot(7,2,k); k=k+1;
stem(fax,fftshift(abs(fft(f7))),'filled')
ftick = get(gca,'XTick');
xlabel('Frequency')
k = k-1;
%set(ax(2:2:end),'XLim',[0 nyq],'XTickLabel',[]);
% set(ax,'XTickLabel',[],'Box','off','TickDir','out');
%set(ax(k),'XTick',[0:10:nyq],'XTickLabel',num2str([0:10:nyq]'))
%set(ax(k-1),'XTick',xtick,'XTickLabel',xticklabel);
%set(ax([1 3 9]),'YLim',[0 1.1]);
%set(ax([5 7]),'YLim',[-1.1 1.1]);
%set(hg,'LineWidth',2);
%set(ax,'YTickLabel',[]);

% Summary of the pairs. 
%
% (1) A gaussian in time is composed of frequencies that fall off from 0 in 
% the shape of a gaussian. This is cool! It's also a unique property of the 
% Gaussian!  Remember that a Fourier transform is just a linear operation. 
% The Gaussian has the special property that it has the same form (Gaussian) 
% under any linear transformation. But the width will be different. If the 
% Gaussian is narrow in time, it will be wide in frequency. And if it is narrow 
% in frequency, it will be broad in time. Guess what: this is related to the uncertainty
% principle! If something is localized in time, it is spread out in
% spectrum (or energy). And vice versa.
%
% (2) A pulse has FT that is periodic under an envelope that
% falls off with frequency. This turns out to be the sin(w)/w function,
% called sinc(w). This is very important as data that is windowed with
% a square envelope is going to show signs of this function. That phenomenon is
% called "ringing".

% (3) A sine function with frequency 1 has an FT with
% weight only at freq=1. We'll deal with phase in a moment. 

% (4) Here's the sine function with frequency 10.

% (5) An exponential looks a lot like an exponential; the weights fall
% off as 1/frequency.

% (6) A delta function is made by adding sinusoids of equal weight at every
% frequency. 

% (7) A comb function has FT that is also a comb function at the
% frequency of the samples and its harmonics.

%
% Now for completeness, let's repeat showing real and imaginary parts.
%
% Run this as a block STARTING HERE...
	figure(5), clf
	k = 1;
	ax(k) = subplot(7,4,k);k=k+1;
	plot(t,f1)
	ax(k) = subplot(7,4,k);k=k+1;
	a = fft(ifftshift(f1));
	plot(fax,fftshift(real(a)));
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(imag(a)))
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(abs(a)))
	set(ax(k-1:-1:k-3),'YLim',1.1*max(abs(a))*[-1 1])
	% stem(fax,fftshift(abs(fft(f1))),'filled')
	ax(k) = subplot(7,4,k);k=k+1;
	plot(t,f2)
	ax(k) = subplot(7,4,k);k=k+1;
	a = fft(ifftshift(f2));
	plot(fax,fftshift(real(a)))
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(imag(a)))
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(abs(a)))
	set(ax(k-1:-1:k-3),'YLim',1.1*max(abs(a))*[-1 1])
	% stem(fax,fftshift(abs(fft(f2))),'filled')
	ax(k) = subplot(7,4,k);k=k+1;
	plot(t,f3)
	ax(k) = subplot(7,4,k);k=k+1;
	a = fft(ifftshift(f3));
	plot(fax,fftshift(real(a)))
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(imag(a)))
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(abs(a)))
	set(ax(k-1:-1:k-3),'YLim',1.1*max(abs(a))*[-1 1])
	ax(k) = subplot(7,4,k);k=k+1;
	plot(t,f4)
	ax(k) = subplot(7,4,k);k=k+1;
	a = fft(ifftshift(f4));
	plot(fax,fftshift(real(a)))
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(imag(a)))
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(abs(a)))
	set(ax(k-1:-1:k-3),'YLim',1.1*max(abs(a))*[-1 1])
	
	ax(k) = subplot(7,4,k);k=k+1;
	plot(t,f5)
	ax(k) = subplot(7,4,k);k=k+1;
	a = fft(ifftshift(f5));
	plot(fax,fftshift(real(a)))
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(imag(a)))
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(abs(a)))
	set(ax(k-1:-1:k-3),'YLim',1.1*max(abs(a))*[-1 1])
	
	ax(k) = subplot(7,4,k);k=k+1;
	plot(t,f6)
	ax(k) = subplot(7,4,k);k=k+1;
	a = fft(ifftshift(f6));
	plot(fax,fftshift(real(a)))
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(imag(a)))
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(abs(a)))
	set(ax(k-1:-1:k-3),'YLim',1.1*max(abs(a))*[-1 1])
	
	ax(k) = subplot(7,4,k);k=k+1;
	plot(t,f7)
	xtick = get(gca,'XTick');
	xlabel('Time [s]')
	ax(k) = subplot(7,4,k);k=k+1;
	a = fft(ifftshift(f7));
	plot(fax,fftshift(real(a)))
	ftick = get(gca,'XTick');
	xlabel('Frequency [Hz]')
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(imag(a)))
	xlabel('Frequency [Hz]')
	ax(k) = subplot(7,4,k);k=k+1;
	plot(fax,fftshift(abs(a)))
	xlabel('Frequency [Hz]')
	set(ax(k-1:-1:k-3),'YLim',1.1*max(abs(a))*[-1 1])
	
	k = k-1;
	set(ax(2:4:end),'XLim',[-nyq nyq]);
	set(ax(3:4:end),'XLim',[-nyq nyq]);
	set(ax(4:4:end),'XLim',[-nyq nyq]);
	set(ax,'Box','off','TickDir','out','Xtick',[])
	set(ax(k:-1:k-2),'XTick',[-nyq 0 nyq])
	set(ax(k:-1:k-2),'XTickLabel',num2str([-nyq 0 nyq]'))
	set(ax(1:3:k-3),'XTick',[0 .5 1],'XTickLabel',[],'TickDir','out')
	set(ax(k-3),'XTick',[0 .5 1],'XTickLabel',num2str([-.5;0;.5]))
	
	set(gcf,'CurrentAxes',ax(1))
	title('s(t)')
	set(gcf,'CurrentAxes',ax(2))
	title('Real part of FT')
	set(gcf,'CurrentAxes',ax(3))
	title('Imaginary part of FT')
	set(gcf,'CurrentAxes',ax(4))
	title('Amplitude of FT')
% ...ENDING HERE

% Pause and reflect. There are a few things to mention. 
%
% First, notice that the zero frequency is special. It represents the mean 
% level: the nonsinusoidal component: the cos(0) component. This component is 
% always purely real.

% Invertibility and Hermitian symmetry. The FT decomposes our signals into sums of
% sines and cosines. This process is invertible. We could start with the
% Fourier transform -- a list of "intensities" at each frequency -- and
% convert it to a list of "intensities" as a function of time. The reason
% we can get away with this is that when we listed our intensities in the
% "frequency domain," we were careful to give two values: amplitude and
% phase, or cosine and sine components, or real and imaginary parts.  To
% ensure that the process is really invertible, we also need to make sure
% that we have the right number of values in the transforms. Put simply, we
% had better represent enough frequencies in the FT so that we can
% reconstitute the original signal by adding together all the right
% sinusoids. Let's look at how many frequencies we need to represent the 
% time functions we've been dealing with. If you left dt=
% 0.01 alone, the answer is 100. 

length(t)
%
% So all of the functions of time that we've considered are described by a
% list this long containing the weights at each time point. What about the
% FT? The highest frequency is the Nyquist, which
% is half the sampling rate. That gives us just 50 frequencies.  How do
% we transform 100 values into just 50 frequencies? You know the answer: at each
% frequency we need 2 numbers, a cosine and a sine component (or amplitude
% and phase). Moreover, for the 0 frequency, there is just one component
% becuase sin(0)=0. Since cos(0)=1, it's clear that the weight at freq=0 is
% just a constant function of time. Also, at the nyquist frequency itself,
% all we can have is a cosine component. To see why, recognize that we only
% have enough samples in time to represent alternating positive and
% negative values at adjacent time points. If you put a 0 at t=0, which you
% must for sin(t), then there's nothing to alternate. Okay, cool, we've got
% two components for each of the frequencies 1 to nyquist-1, one weight at
% freq=0 and one weight at freq=nyquist. That's 100 values in the FT, which
% corresponds to 100 values in the time function. Perfect!
%
% So what's the deal with the negative frequencies? The
% frequency axis does not go from 0 to nyquist but from -nyquist to
% nyquist-1. Doesn't that screw up our bookkeeping?  Here's the answer. We
% have only considered the Fourier transform of real valued functions of time
% (left column of Fig. 5). The Fourier transform returns complex numbers at
% each frequency, which we can think of a sine waves (imaginary part) and
% cosine waves (real part). In fact the Fourier transform can act
% on complex functions of time. Instead of a real value at each time, we
% could have a complex number. We may not have any use for this, but the
% full transform takes these 2-vectors (real and imaginary parts) and
% transforms them to 4 values: a real and imaginary value at two
% frequencies, one positive and one negative. Of course, there are half as
% many frequencies as there are time points. What I'm saying is that we
% actually had twice as many values in our time functions as we thought. We
% used real functions of time. The extra 100 values were zeros: the weights
% on the imaginary part of the signals.  So the transform is 1 to 1.
%
% That's all fine and good in theory, but what we care about for most
% applications is the Fourier transform of real valued functions. That's
% like saying that we know that the imaginary part of the functions in the
% left column of Fig. 5 is 0 at every time point. So, I'll ask again, why
% should the FT have 200 values to represent a signal that we can describe
% by 100 points in time? The answer is that the FTs of real signals really
% can be represented with just 100 numbers, just like we said. You really
% don't need both positive and negative frequencies. If I tell you the
% weights for one of them (e.g., +5 Hz component), you know the weights of
% the negative frequency (e.g., -5 Hz component). This shows up as a kind
% of symmetry in the graphs in Figure 5. See if you can deduce the rule from 
% these examples.
%
% I will pose this as a homework problem. First a reminder. An even symmetric
% function obeys the rule f(x) = f(-x). An example is cos(x). An odd symmetric
% function is one that obeys the rule f(x) = -f(-x). An example is sin(x).  
%
% QUESTION: Which of the seven functions in Fig. 5 are even? Which is odd?
% Which is neither even nor odd? What do you notice about the Fourier
% transforms of the even and odd symmetric functions. Answer this for the
% three categories: even, odd, neither even nor odd.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part II. Convolution
%
% One of the most important uses of the FT is for computing convolution.
% Convolution is an operation that looks like filtering. It takes a signal
% and converts it to a new signal. Usually, we can think of the new signal
% as a filtered version of the old one. We gain a deep intuition for the
% process by conceptualizing the operation in the time domain and in the
% frequency domain (i.e., via the fourier transform). I think this is one
% of the most beautiful concepts in applied math -- and it comes up all the
% time (or often if you prefer a frequency based description!).
%
% **** II.1 Convolution in the time domain. ****
%
% In this section, we spend a lot of time understanding convolution without
% talking about FTs. But I want you to keep FTs in the back of your mind
% while we do this. What you should be asking yourself is this. When
% filtering looks like blurring, isn't that a lot like attenuating sharp
% aspects of a signal, and isn't that a lot like attenuating high
% frequencies?
%
% Let's get started.
% Consider the simplest of all possible signals: a delta function,
% otherwise known as a click (if it were sound pressure as a function of
% time) or a line (if it were light intensity as a function of horizontal
% position).

timeOfImpulse = .1
s1 = zeros(size(t)); s1(t==timeOfImpulse)=1;
figure(1),clf,plot(t,s1)

% Consider a blurring function. Let's start with one that you've grown
% accustomed to: eponential decay. We'll see in a moment that this function
% is causal, a useful property when it comes to time-dependent functions.
tau = .03;      % choose a time constant
b = exp(-t/tau);
plot(t,b)
set(gca,'XLim',[0 timeOfImpulse+10*tau])
xlabel('Time [out to 10\tau]')

% Please run this next step without thinking about it. It's a trick that
% allows me to make some pictures easily. You can come back to this if you
% want to pursue the matrix version of convolution.
A = tril(toeplitz(b));
% figure(2), clf, imagesc(A)
% plot(A(:,end))

% What is convolution?
%
% It is a blurring of s1 by b. It is a new function of time that can be
% thought of as a weighted sum of s1. It's actually easy to see what it is
% by looking at the formula.  
%
% s2(k) = Sum_over_d{product of s1(d) * b(k-d)}
%
% s2 is the result of the convolution. It is a function of time, but you've
% probably noticed that the variable t is not in the equation.  Actually
% both k and d refer to time, but because we are about to shift things
% around time takes on several meanings.   It's actually easier to keep
% things straight if we keep t out of the equation. The convolution is
% defined at each time, which I will denote this by k. You will come to
% view k as an offset of the blur function.  The original signal remains a
% function of time, but we'll use a dummy variable for time. That's what d
% stands for: we can refer to d as "dime".  There is nothing special about
% s1(d). It is the signal plotted as a function of "dime" instead of time.
% Same thing -- no problem. What is b(k-d)?  It is b reflected about the
% y-axis and then offset by k time steps to the right.  Here is what b(k-d)
% looks like when k=0.
d = t;
figure(1), clf
plot(-d,b)
xlabel('-dime')
% Here is what b(k-d) looks like for increasing values of k. 
figure(1),clf
ax = [];
for k = 1:10
    ax(k) = subplot(10,1,k);
    plot(t(k)-d,b)
    set(gca,'XLim',[-10*tau timeOfImpulse+10*tau],'Box','off')
end
xtick = get(ax(1),'XTickLabel');
ytick = get(ax(1),'YTickLabel');
set(ax,'XTickLabel',[],'YTickLabel',[]);
set(ax(k),'XTickLabel',xtick,'YTickLabel',ytick)
xlabel('Time [d]')
set(gcf,'CurrentAxes',ax(1))
title('b(k-d) for k=1 to 10') 
% It is a picture of the function reflected around the ordinate and then
% offet to the right by k. Can you see this? 

% Now let's generate the convolution, s2(k), step of k by step of k.
% First look at the equation and convince yourself that there is one value of the
% convolution produced at each offset step. Run the animation in the next
% loop. You might want to put some pauses in the code. Here's what you
% should see. In the top panel, s1(d) is shown by a stem function in red. It is
% just an impulse. The blur function, b(k-d) is shown in blue. It is
% stepped along as k increases. The middle panel shows the product of s1(d)
% times b(k-d). Notice that there is a value at each value for dummy time.
% The bottom plot shows the sum of the products. At each time step, k, we
% get just one value. This is the convolution, s2(k). You will notice that
% the sum is weighted by the width of the time bins. Why is this? What
% effect does this have on the relationship between s1 and s2?
figure(2),clf
h = [];
h(1) = subplot(3,1,1); hold on;
h(2) = subplot(3,1,2);
h(3) = subplot(3,1,3); hold on;
% set(h,'TickDir','out','Box','off')
j = min([min(find(t>.1 + 10*tau)) length(t)-1])
for k = 1:j
    set(gcf,'CurrentAxes',h(1))
    if k>1, delete(hp), end
    stem(d,s1,'r'), hold on;
    hp = plot(d,flipud(A(:,end-k+1)));
    hold off
    
    set(gcf,'CurrentAxes', h(2))
    s1TimesConvKernel = s1 .* flipud(A(:,end-k+1));
    stem(d,s1TimesConvKernel,'filled')
   
    set(gcf,'CurrentAxes', h(3))
    stem(t(k),sum(s1TimesConvKernel));
    set(h,'XLim',[-.05 timeOfImpulse+10*tau])
    set(h(2:3),'YLim',[0 1])
    % set(h(3),'YLim',[0 .02])
    set(h,'Box','off','TickDir','out')
    drawnow
end

% This is the original function, a delta function or impulse, smeared out
% by the exponential blurring function. The technical term for smearing is
% filtering! What you've simulated is a very fine click played through a
% woofer or a thin voltage spike as seen through an RC filter.  

% Oh yeah, there's an easy way to get matlab to do the convolution for you.
% Type 'help conv' in the command window and read what it says. Here's our
% friend, the impulse convolved with the exponential. Notice that the
% result is vector whose length is bigger than t.  
s2 = conv(s1, b);
figure(3),clf
hold on
stem(t,s1,'r'), plot(t,b,'g',t,s2(1:length(t)),'r--')
% Here's the convolution between a comb function and a guassian
s2 = conv(comb10,gaus1);
clf, hold on
stem(t,comb10,'k'),plot(t,gaus1,'k',t,s2(1:length(t)),'r--')

% Try playing with other functions

% **** II.2 Convolution and Fourier transforms ****
% 
% Fourier transforms make convolution easy. That's because convolution
% between two functions of time, like s1(t) and b(t), is equivalent to
% multiplying their fourier transforms. To be more precise, if S1(w) is the
% fourier transform of s1(t) and B(w) is the FT of b(t), then
%
% S2(w) = S1(w) .* B(w) 
%
% is the FT of s2(t).
%
% This provides a simple recipe for computing convolutions. (1) Take the
% fourier transforms of the two functions you want to convolve. (2)
% multiply the FTs (there's a complex number at each frequency!). (3) Take
% the inverse fourier transform. 
%
% Let's do this with the exponential blurring function and s1. To work with
% b, I need it to be defined on the same time axis as s1. I'll call it
% newb.
newb = zeros(size(t));
newb(1:length(b)) = b;
% Take the fourier transforms
B = fft(newb);
S1 = fft(s1);
% Multiply them at each frequency.
S2 = B .* S1;
% Take the inverse fourier transform
s2 = ifft(S2);

% How does s2 look?
% You have already seen the fourier transforms of the signal and
% exponential. So I'm just going to plot the result
figure(4), clf
plot(t,s2)

% Let's illustrate with a more realistic set of functions.
s1 = sin(2*pi*3*t) + .33 * sin(2*pi*9*t) + .2 * sin(2.*pi*15*t);
% Corrupt this signal by adding noise
s1 = s1 + rand(size(s1));
S1 = fft(s1);
% You could use a different blurring function. I'll apply a gaussian blur.
% See if you can figure out what the fftshift is doing here.
g = fftshift(exp(-((t-.5)/tau).^2));
G = fft(g);
S2 = G .* S1;
s2 = ifft(S2);
figure(4), clf
subplot(3,2,1)
plot(t,s1)
subplot(3,2,2)
plot(fax, fftshift(abs(S1)));
subplot(3,2,3)
plot(t,g)
subplot(3,2,4)
plot(fax,fftshift(abs(G)))
subplot(3,2,5)
plot(t,s2)
subplot(3,2,6)
plot(fax, fftshift(abs(S2)));

return

% Additional topics. You should work on all of these. If they paid me more to
% teach this class, I would continue the tutorial to cover these topics.
% But I'll leave them to you.
%
% 1. Develop better intutions about filters and multiplication in the freq domain.
% Find real world examples: muffling of a cheap speaker, resonance and timbre in
% musical instruments, blurring of an image by an imperfect lense. 
%
% 2. What is a bode plot? What is a power spectrum?
%
% 3. What is an impulse response function? Why is it useful? If I play a
% click through a speaker, it comes out a little less sharp. The shape of
% this output function is an example of an impulse response function. What
% is it good for? How does it relate to what your speaker will do when you
% send it signals that represent music? Answer this question in terms of
% convolution and fourier transforms. What are the critical assumptions in
% this argument (hint, there are 2 tenets to linear systems theory). 
% 
% 4. Filtering a signal is equivalent to attenuating frequencies. Given
% this idea, think about designing a filter by specifying its amplitude
% spectrum. What problems do you face? Why can't you make a filter that
% passes a band of frequencies from 10-20 Hz but cuts out all other
% frequencies? (Hint1: FT pairs are symmetrical if you know f(t) <---> F(w)
% are a transform pair, then F(t) <---> f(w) are also transform pairs. Hint
% 2: apply hint 1 to the pulse(t)<-->sinc(w) transform pair.)
%
% 5. Sampling & aliasing. Consider a signal with high frequencies in it. 
% y = sin(2*pi*3*t) + .33 * cos(2*pi*9*t) + .2 * sin(2.*pi*15*t);  
y = sign(cos(2*pi*27*t - pi/7));
figure(1), plot(t,y)
% Suppose the signal is sampled not at every time point but only once every
% 100 msec
ysampled = y .* comb10;
plot(t,ysampled)
% Why does the sampled waveform look so crumby? Look at the FT of ysampled.
% Notice the amplitude at low frequencies. What are they doing there? We
% did not make y using low frequencies? They are "aliases". Recall that
% convolution between two time-functions corresponds to multiplication
% between their FTs. What do you think multiplying two time functions
% corresponds to?  With the answer to that question in mind, what is the
% fourier transform of y, comb10 and what do you think you would do with
% these fourier transforms to produce the fourier transform of ysampled?
% Think hard about what is going on in these graphs.  Why do you think this
% undersampling of y produces "aliasing"? What is going on in the frequency
% domain when you multiply y by comb10?











