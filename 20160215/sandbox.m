%% construct a 10Hz signal, sampled at 1000Hz. 
% Recalling that the frequency f of a sine wave is given by y=sin(2pift)
% Notice the general approach in defining time series data: we first 
% construct a timebase (conventionally named tvec, or t) and then the signal
%  first construct a 10Hz signal, sampled at 1000Hz

Fs1 = 1000; % Fs is the conventional variable name for sampling freq
F1 = 10; twin = [0 1]; % use a 1-second time window (from 0 to 1s)
tvec1 = twin(1):1/Fs1:twin(2); % timebase for signal
signal1 = sin(2*pi*F1*tvec1);
%% plot the 10hz signal
plot(tvec1,signal1);

%% Let's say we are going to sample this signal at 12Hz:
Fs2 = 12;
tvec2 = twin(1):1/Fs2:twin(2);
% we have an input signal, specified by tvec1 and signal1; 
% return those values of signal1 for those values of tvec1 closest to 
% those in tvec2.
signal2 = interp1(tvec1,signal1,tvec2,'nearest');

Fs3 = 20;
tvec3 = tvec1(1:1000/Fs3:end);
signal3 = signal1(1:1000/Fs3:end);

% plot signal 1 and signal 2
plot(tvec1,signal1);
hold on;
plot(tvec2,signal2,'.g','MarkerSize',20); 
plot(tvec1,-sin(2*pi*2*tvec1),'r--','LineWidth',2); % 2hz signal:aliasing (i.e. the existence of multiple underlying signals which could produce a given set of samples), 
%stem(tvec3,signal3,'.g','MarkerSize',20);
xlabel('time (s)'); ylabel('y');

%% Subsampling (decimating) time series data
% generate a signal consisting of two frequencies
Fs1 = 1200;
F1 = 3; F2 = 10;
twin = [0 1];
 
tvec1 = twin(1):1/Fs1:twin(2);
signal1a = sin(2*pi*F1*tvec1); signal1b = 0.5*sin(2*pi*F2*tvec1);
signal1 = signal1a + signal1b;

dt = 100;
tvec2 = tvec1(1:dt:end);
signal2 = signal1(1:dt:end); % sample at 12 Hz - every 100th sample
 
subplot(131)
plot(tvec1,signal1);
hold on;
stem(tvec2,signal2,'.g','MarkerSize',20);
title('without anti-aliasing filter');

%% sample at 12 Hz with different method
% The way decimate() works is that it first applies a filter to the data, 
% removing any frequencies that could cause aliases (i.e. anything with a 
% frequency of at least half the new sampling frequency). 
tvec1d = decimate(tvec1, dt);
signal2d = decimate(signal1,dt);
 
subplot(132)
plot(tvec1,signal1a,'b--');
hold on;
stem(tvec1d,signal2d,'.g','MarkerSize',20);
xlabel('time (s)'); ylabel('y');
title('with anti-aliasing filter');
 
subplot(133)
plot(tvec2,signal2-signal2d,'r','LineWidth',2);
title('difference');

%% Reconstructing a signal from sampled data (optional)
% constructing a 100 Hz signal, sampled at 2 kHz
fs1 = 2000;
tvec = 0:1/fs1:4; % construct time axis, sampled at 2kHz
 
freq1 = 100;
y = sin(2*pi*freq1*tvec); % 100Hz signal
 
ax1 = subplot(211);
stem(tvec,y); title('original');

% Say we wish to subsample down to 500Hz to save space. 
% Naively, we might simply take every 4th sample (in fact, 
% this is what the MATLAB function downsample() does):
subsample_factor = 4;
 
tvec2 = tvec(1:subsample_factor:end); % take every 4th sample
y2 = y(1:subsample_factor:end);
 
ax2 = subplot(212);
stem(tvec2,y2,'r'); title('subsampled');
xlabel('time (s)');

xl = [1 1.04];
linkaxes([ax1, ax2], 'x');
set(ax1,'XLim',xl); % see what I did there?)


hold on;
 
y_interp = interp1(tvec2,y2,tvec,'linear');
p1 = plot(tvec,y_interp,'b');
 
y_interp2 = interp1(tvec2,y2,tvec,'spline');
p2 = plot(tvec,y_interp2,'g');
 
legend([p1 p2],{'linear','spline'},'Location','Northeast'); legend boxoff





