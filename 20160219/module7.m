%% basics using spectrogram() function
% cd to your location here
rootDir = '/Users/sirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R016/R016-2012-10-03';
cd(rootDir);
cfg = [];
cfg.fc = {'R016-2012-10-03-CSC04a.Ncs'};
csc = LoadCSC(cfg);
Fs = csc.cfg.hdr{1}.SamplingFrequency;
%% restrict data
cscR = restrict(csc,3282,3286); % if you don't have this, implement it (it's one line of code!)
plot(cscR.tvec,cscR.data); % note there are various gamma oscillations present, as well as a large negative-going transient
%% construct and plot the spectrogram
wSize = 512;
WINDOW = hanning(wSize);
NFFT = 1:200;
NOVERLAP = wSize/2; % step size = wSize - NOVERLAP;
% stepsPerSec = round(Fs/(wSize - NOVERLAP);
[S,F,T,P] = spectrogram(cscR.data,WINDOW,NOVERLAP,NFFT,Fs); % [S,F,T,P] = spectrogram(X,WINDOW,NOVERLAP,NFFT,Fs)
figure
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
%% raw LFP overlaid onto the spectrogram
hold on;
lfp_minmax = 25; lfp_cent = 125; % range and mean of LFP plotting
tvec0 = cscR.tvec - cscR.tvec(1); % align LFP with spectrogram make cscR.tvec(1) to be zero
data = rescale(cscR.data,-lfp_minmax,lfp_minmax); % rescale so easier to visualize
data = data+lfp_cent;
lfp_h = plot(tvec0,data,'k');

%% spectrogram parameters
% change several parameters that do not affect the power estimate
wSize = 512;
WINDOW = hanning(wSize);
NFFT = 1:.25:200;
%NFFT = 1:.01:200;
NOVERLAP = 384; % step size = wSize - NOVERLAP;
%NOVERLAP = wSize - 1; % ultra-smooth spectrogram
% stepsPerSec = round(Fs/(wSize - NOVERLAP);
[S,F,T,P] = spectrogram(cscR.data,WINDOW,NOVERLAP,NFFT,Fs); % [S,F,T,P] = spectrogram(X,WINDOW,NOVERLAP,NFFT,Fs)
figure
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  

%% change shape of the window
wSize = 512;
WINDOW = rectwin(wSize);
NFFT = 1:.25:200;
%NFFT = 1:.01:200;
NOVERLAP = 384; % step size = wSize - NOVERLAP;
%NOVERLAP = wSize - 1; % ultra-smooth spectrogram
% stepsPerSec = round(Fs/(wSize - NOVERLAP);
[S,F,T,P] = spectrogram(cscR.data,WINDOW,NOVERLAP,NFFT,Fs); % [S,F,T,P] = spectrogram(X,WINDOW,NOVERLAP,NFFT,Fs)
figure
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  