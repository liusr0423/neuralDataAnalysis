%% construct a 10Hz signal, sampled at 1000Hz. 
% Recalling that the frequency f of a sine wave is given by y=sin(2?ft)
% Notice the general approach in defining time series data: we first 
% construct a timebase (conventionally named tvec, or t) and then the signal

Fs1 = 1000; % Fs is the conventional variable name for sampling freq
F1 = 10; twin = [0 1]; % use a 1-second time window (from 0 to 1s)
tvec1 = twin(1):1/Fs1:twin(2); % timebase for signal
signal1 = sin(2*pi*F1*tvec1);

% Let's say we are going to sample this signal at 12Hz:
Fs2 = 12;
tvec2 = twin(1):1/Fs2:twin(2);
% we have an input signal, specified by tvec1 and signal1; 
% return those values of signal1 for those values of tvec1 closest to 
% those in tvec2.
signal2 = interp1(tvec1,signal1,tvec2,'nearest');

% plot signal 1 and signal 2
plot(tvec1,signal1);
hold on;
plot(tvec2,signal2,'.g','MarkerSize',20);
plot(tvec1,-sin(2*pi*2*tvec1),'r--','LineWidth',2);
xlabel('time (s)'); ylabel('y');