%% plot a simple sinusoid
Fs = 100; % in samples per second (Hz)
t0 = 0; t1 = 1; % start and end times
tvec = t0:1./Fs:t1; % construct time axis
 
f = 2; % frequency of sine to plot: 2Hz
y = sin(2*pi*f*tvec); % note sin() expects arguments in radians, not degrees (see also ''sind()'')
 
stem(tvec,y);

%%  sinusoid is further characterized by its phase and amplitude
phi = pi/2;
 
subplot(221)
y = sin(2*pi*f*tvec + phi); % a phase shift
stem(tvec,y);

hold on;
plot(tvec,cos(2*pi*f*tvec),'r--','LineWidth',2) % notice, cosine is simply phase shifted sine
legend('sin (phase-shifted)', 'cos');
 
a = 2;
subplot(222)
y = a.*sin(2*pi*f*tvec + phi); % amplitude change
stem(tvec,y); % note scale of y axis!

%% frequency modulation (FM) signal
% cross-frequency coupling
f2 = 10;
m = 2;
 
subplot(311)
s1 = sin(2*pi*f*tvec);
plot(tvec,s1); title('message');
 
subplot(312);
s2 = sin(2*pi*f2*tvec);
plot(tvec,s2); title('carrier');
 
subplot(313);
% the phase of one (slower, ?message?) signal is related (?coupled?) to 
% the frequency of another (faster, ?carrier?) signal.
s3 = sin(2*pi*f2*tvec + m.*sin(2*pi*f*tvec - pi/2));
plot(tvec,s3); title('FM signal');

%% harmonic series example
mag = [0.1 0 1.3 0.5]; % magnitudes for each term
pha = [-pi/6 0 pi 2*pi/3]; % phases for each term
f = 2; % base frequency
signal_out = zeros(size(tvec));
for ii = 1:numel(mag) % note, the book chapter uses i, not best practice!
 
    this_signal = mag(ii)*cos(2*pi*f*ii*tvec + pha(ii));
    plot(tvec,this_signal,'r:'); hold on;
    signal_out = signal_out + this_signal; % build the sum
 
end
plot(tvec,signal_out,'LineWidth',2);

% In fact, the central insight underlying Fourier analysis is that we can 
% use a sum (series) of sinusoids to approximate any signal to arbitrary 
% precision. An important corollary to this is that any signal can be 
% decomposed into a series of sinusoids


%% Decomposing and reconstructing a signal
rng('default'); % reset random number generator to reproducible state, so your plot will look like mine!
x = round(rand(1,8)*10); % generate a length 8 vector of integers between 0 and 10
xlen = length(x);

% get magnitudes and phases of Fourier series
X = fft(x);
Xmag = abs(X); % magnitudes, a_n
Xphase = angle(X); % phases, phi_n

n = 0:xlen-1;
t = 0:0.05:xlen-1; % a finer timescale to show the smooth signal later
 
for iH = xlen-1:-1:0 % reconstruct each harmonic
    s(iH+1,:) = Xmag(iH+1)*cos(2*pi*n*iH/xlen + Xphase(iH+1))/xlen;
    sm(iH+1,:) = Xmag(iH+1)*cos(2*pi*t*iH/xlen + Xphase(iH+1))/xlen;
    % detail: xlen appears here because the fundamental frequency used by fft() depends on this
end

ssum = sum(s); % coarse timescale (original points)
smsum = sum(sm); % fine timescale (to see full signal)
 
figure;
plot(n, x, 'go', t, smsum, 'b', n, ssum, 'r*');
legend({'original','sum - all','sum - points only'});

%% Interpreting the output of MATLAB's fft() function
Fs = 20; % in samples per second (Hz)
t0 = 0; t1 = 1; % start and end times
tvec = t0:1/Fs:t1-(1/Fs); % construct time axis; generate exactly 20 samples
 
f = 2; % signal frequency
y = sin(2*pi*f*tvec); % construct signal, a 2Hz sine wave sampled at 20Hz for 1s
 
yfft = fft(y,length(y));
yfft_mag = abs(yfft); yfft_ph = angle(yfft);
stem(yfft_mag)

%% modify using fft()
Npoints = length(y);
F = [-Npoints/2:Npoints/2-1]./Npoints; % construct frequency axis
 
yfft_mag = fftshift(yfft_mag); % align output,fftshift() cuts the second 
% (complex) half of the spectrum and pastes it backwards at the beginning, 
% so that our frequency axis is now correct;it is in units of 1 / Fs so 0.1 
% corresponds to the 2Hz we put in.

stem(F,yfft_mag);
 
xlabel('Frequency (Fs^{-1})');

%% zero-padding
tvec = t0:1/Fs:t1; % Change the tvec variable above to contain one more 
                   % sample so our signal is now no longer an integer 
                   % number of periods
nPoints = [length(tvec) 64 256 1024];

for iP = 1:length(nPoints) % repeat fft with different numbers of points
 
    nP = nPoints(iP);
    subplot(2,2,iP);
  
    y = sin(2*pi*f*tvec);
    yfft = fft(y,nP); %using the second, optional, argument of fft() to zero pad the signal
    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
 
    F = [-nP/2:nP/2-1]./nP;
    yfft_mag = fftshift(yfft_mag);
    plot(F,yfft_mag,'kx',F,yfft_mag,'k');
 
    title(sprintf('%d point FFT',nP));
    xlabel('Frequency (Fs^{-1})');
 
end
% A typical value to use for the number of points to evaluate the FFT is 
% the next power of 2 (after however many samples your signal contains). 
% This is because part of what makes the FFT fast is that it can easily 
% divide the signal in half. But non-power of 2 values also work.

%  zero-padding the original 21-point signal rather than making the signal
%  longer
tvec = t0:1/Fs:t1; % our signal is now no longer an integer number of periods
y =  [sin(2*pi*f*tvec),zeros(1,1024-length(tvec))]; % pad signal with zeros to 1024 points 
nP = length(y);
yfft = fft(y,nP);
yfft_mag = abs(yfft); yfft_ph = angle(yfft); 
F = [-nP/2:nP/2-1]./nP;
yfft_mag = fftshift(yfft_mag);
plot(F,yfft_mag,'kx',F,yfft_mag,'k');


%% Spectral leakage
% increase the length of the signal by repeating it a number of times
tvec = t0:1/Fs:t1;
nRepeats = [1 2 4 8];
 
nP =  1024; % use 1024 points to evaluate the FFT 
 
for iP = 1:length(nRepeats)
 
    subplot(2,2,iP);
 
    y = sin(2*pi*f*tvec);
    y = repmat(y,[1 nRepeats(iP)]); % repeat the signal a number of times
 
    yfft = fft(y,nP); % zero-pad to nP
    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
 
    F = [-nP/2:nP/2-1]./nP;
    yfft_mag = fftshift(yfft_mag);
    plot(F,yfft_mag,'kx',F,yfft_mag,'k');
 
    title(sprintf('%d repeats',nRepeats(iP)));
    xlabel('Frequency (Fs^{-1})');
 
end

%% windowing
nP = 25;
nPFFT = 1024;
 
windows = {'rectwin','triang','hamming','hanning','blackman'};
cols = 'rgbcmyk';
 
for iW = 1:length(windows)
 
    eval(cat(2,'wn = ',windows{iW},'(nP);')); % make sure you understand this
    wn = wn./sum(wn); % normalization to make sure integral of the window is 1 to preserve power estimates
 
    subplot(211); % plot the window
    plot(wn,cols(iW),'LineWidth',2); hold on;
 
    subplot(212);
    yfft = fft(wn,nPFFT);
    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
 
    F = [-nPFFT/2:nPFFT/2-1]./nPFFT;
    yfft_mag = fftshift(yfft_mag);
 
    h(iW) = plot(F,yfft_mag,cols(iW),'LineWidth',2); hold on;
 
end
 
xlabel('Frequency (Fs^{-1})');
legend(h,windows);
set(gca,'XScale','log');
%%
% the spectrum of the windowed signal equals the convolution of the 
% signal's spectrum and the window's spectrum;
% multiplication in the time domain equals convolution in the frequency
% domain;

% use a Hamming window instead of a rectangular window on spectral leakage
tvec = t0:1/Fs:t1;
%nRepeats = 4 ;
 
nP =  1024; % use 1024 points to evaluate the FFT 
 
y = sin(2*pi*f*tvec);
wn = hamming(length(tvec))';
wn = wn ./sum(wn);
y2 = y .* wn;  % Hamming window
%y = repmat(y,[1 nRepeats]); % rectangular window: repeat the signal a number of times

% subplot(2,1,1);
% yfft = fft(y,nP); % zero-pad to nP
% yfft_mag = abs(yfft); yfft_ph = angle(yfft);
% F = [-nP/2:nP/2-1]./nP;
% yfft_mag = fftshift(yfft_mag);
% plot(F,yfft_mag,'kx',F,yfft_mag,'k');
% title(sprintf('%d repeats',nRepeats));
% xlabel('Frequency (Fs^{-1})');
    
%subplot(2,1,2);
yfft = fft(y2,nP); % zero-pad to nP
yfft_mag = abs(yfft); yfft_ph = angle(yfft);
F = [-nP/2:nP/2-1]./nP;
yfft_mag = fftshift(yfft_mag);
plot(F,yfft_mag,'kx',F,yfft_mag,'k');
title(sprintf('Hamming window'));
xlabel('Frequency (Fs^{-1})');




%% Robust spectral estimation methods
% Power Spectral Density (PSD): Pxx = periodogram(X,WINDOW)
[Pxx,F] = periodogram(y,[],nP,Fs); % X is windowed with a rectangular window by default
plot(F,Pxx); xlabel('Frequency (Hz)');
ylabel('Power');
hold on;
[Pxx,F] = periodogram(y,hanning(length(y)),nP,Fs);
plot(F,Pxx,'r');

% this method is biased: its variance does not tend to zero as data length
% goes to infinity; this makes spectrum look noisy
%% Welch's method
% cut the signal up into smaller segments, estimate the spectrum for 
% each segment, and combine the estimates (Bartlett's method)
% Welch's method uses segments (or windows) that can overlap.

Fs = 20; % in samples per second (Hz)
t0 = 0; t1 = 1;
f = 2;
nRepeats = 4;
 
tvec = t0:1/Fs:t1-(1/Fs);
 
nP =  1024;
y = sin(2*pi*f*tvec);
y = repmat(y,[1 nRepeats]);
 
[Pxx,F] = periodogram(y,rectwin(length(y)),nP,Fs);
plot(F,Pxx);
 
hold on;
wSize = 40;
%  [Pxx,F] = pwelch(X,WINDOW,NOVERLAP,F,Fs)
%  uses NOVERLAP samples of overlap from section to section
[Pxx,F] = pwelch(y,rectwin(wSize),wSize/2,nP,Fs);
plot(F,Pxx,'r'); xlabel('Frequency (Hz)');



%% Pitfalls for real-world signals
% Neuralynx data: chunks of 512 samples are sampled at a regular interval   
% (Fs = 2000) but the interval between chunks is slightly larger than 
% within, such that the overall Fs is slightly smaller than 2000. 


Fs = 20; % in samples per second (Hz)
t0 = 0; t1 = 1; f = 2;
nP =  1024;
gaps = [5 10 15]; % idx of samples to be removed
 
tvec = t0:1/Fs:t1;%-(1/Fs);
y = sin(2*pi*f*tvec);

subplot(211)
plot(tvec,y,'k*'); hold on;

yfft = fft(y,nP);
yfft_mag = abs(yfft); yfft_ph = angle(yfft);
 
F = [-nP/2:nP/2-1]./nP;
yfft_mag = fftshift(yfft_mag);
 
subplot(212);
plot(F,yfft_mag,'k-x'); hold on;
 
xlabel('Frequency (Fs^{-1})');

% signal with gaps
y = sin(2*pi*f*tvec);
y2 = y;
y2(gaps) = []; tvec(gaps) = []; % remove
 
subplot(211);
plot(tvec,y2,'ro'); hold on;
legend('signal','signal(with gas)')

yfft = fft(y2,nP);
yfft_mag = abs(yfft); yfft_ph = angle(yfft);
 
F = [-nP/2:nP/2-1]./nP;
yfft_mag = fftshift(yfft_mag);
 
subplot(212);
plot(F,yfft_mag,'b-x');
legend('y','y2 (with gaps)')

% the estimate for the data with some samples removed (in blue) is 
% substantially lower than the true frequency.

%% Application to real data
% cd to your location here
rootDir = '/Users/sirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R016/R016-2012-10-08_promoted';
cd(rootDir);
% loading a ventral striatal LFP signal
cfg = [];
cfg.fc = {'R016-2012-10-08-CSC04d.ncs'};
csc = LoadCSC(cfg);
%% first rest segment
% restrict to prerecord, leaving some time (10s) before rat actually goes on track
csc_pre = restrict(csc,0,csc.cfg.ExpKeys.TimeOnTrack(1)-10);
%% check if sampling is ok
plot(diff(csc_pre.tvec)); % only minimal differences
Fs = 1./mean(diff(csc_pre.tvec));
%% decimate to speed up processing
dsf = 4;
csc_pre.data = decimate(csc_pre.data,dsf);
csc_pre.tvec = downsample(csc_pre.tvec,dsf);
csc_pre.cfg.hdr{1}.SamplingFrequency = csc_pre.cfg.hdr{1}.SamplingFrequency./dsf;
Fs = 1./mean(diff(csc_pre.tvec));

%% compute the spectrum
[Pxx,F] = periodogram(csc_pre.data,hamming(length(csc_pre.data)),length(csc_pre.data),Fs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);

%% use Welch spectrum instead, using a Hamming window of the same size and 50% overlap
wSize = 1024;
[Pxx,F] = pwelch(csc_pre.data,hamming(wSize),[],[],Fs);
hold on;
plot(F,10*log10(Pxx),'r'); xlabel('Frequency (Hz)');ylabel('Power (dB)');


