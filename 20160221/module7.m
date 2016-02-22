%% basics using spectrogram() function
% cd to your location here
rootDir = '/Users/liusirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R016-2012-10-03';
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
% increasing the frequency and time resolution in this way does not 
% change the underlying estimate since the spectral content is estimated 
% over a window that is 512 samples wide.
%% changing type of the window
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
%% compare spectrogram with different window size
wSizes = [512;1024];
for n = 1:length(wSizes)
    WINDOW = hanning(wSizes(n));
    NFFT = 1:200;
    %NFFT = 1:.01:200;
    NOVERLAP = wSizes(n)/2; % step size = wSize - NOVERLAP;
    %NOVERLAP = wSize - 1; % ultra-smooth spectrogram
    % stepsPerSec = round(Fs/(wSize - NOVERLAP);
    [S,F,T,P] = spectrogram(cscR.data,WINDOW,NOVERLAP,NFFT,Fs); % [S,F,T,P] = spectrogram(X,WINDOW,NOVERLAP,NFFT,Fs)
    subplot(length(wSizes),1,n)
    imagesc(T,F,10*log10(P)); % converting to dB as usual
    set(gca,'FontSize',20);
    axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');
end
% the larger window has tighter frequency estimates, but the spectrogram is 
% now more smeared out in time
%% mind the gaps
cscR = restrict(csc,3300,3340);
 
[S,F,T,P] = spectrogram(cscR.data,rectwin(256),128,1:200,Fs);
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
 
hold on;
 
lfp_minmax = 25; lfp_cent = 125; % range and mean of LFP plotting
tvec0 = cscR.tvec - cscR.tvec(1); % align LFP with spectrogram
data = rescale(cscR.data,-lfp_minmax,lfp_minmax); 
data = data+lfp_cent;
 
lfp_h = plot(tvec0,data,'k');
xlim([tvec0(1) tvec0(end)]);

% spectrogram() ignores the gap present in the LFP, causing the spectrogram 
% and the LFP to become misaligned
%% Processing Neuralynx events files (.nev)
cfg = [];
cfg.eventList = {'Feeder 0 nosepoke','Feeder 1 nosepoke', ...
    '1 pellet cue','3 pellet cue','5 pellet cue', ...
    '1 pellet dispensed','3 pellet dispensed','5 pellet dispensed'};
cfg.eventLabel = {'n0','n1', ...
    'c1','c3','c5', ...
    'd1','d3','d5'}; 
evt = LoadEvents(cfg);

%% Loading Neuralynx data into FieldTrip
% remember to cd to your data folder
rootDir = '/Users/liusirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R016-2012-10-03';
cd(rootDir);
fc = {'R016-2012-10-03-CSC02b.ncs'};
data = ft_read_neuralynx_interp(fc);

cfg = [];
cfg.t = cat(2,getd(evt,'n0'),getd(evt,'n1'));
cfg.mode = 'nlx';
cfg.hdr = data.hdr;
cfg.twin = [-1 4];
 
trl = ft_maketrl(cfg);
 
cfg = [];
cfg.trl = trl;
data_trl = ft_redefinetrial(cfg,data); 

%% Constructing the event-triggered spectrogram
cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.channel      = 'R016-2012-10-03-CSC04a';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 10:2:100; % frequencies of interest: 10hz to 100hz in steps of 2Hz(1/windowSizeinSec)
cfg.t_ftimwin    = ones(size(cfg.foi)).*0.5;  % window size: fixed at 0.5s
cfg.toi          = -0.5:0.05:3.5; % times of interest: time window "slides" 
                                  % from -0.5 to 3.5 sec in steps of 0.05 sec (50 ms)
 
TFR = ft_freqanalysis(cfg, data_trl);
 
figure
cfg = []; cfg.channel = 'R016-2012-10-03-CSC04a';
ft_singleplotTFR(cfg, TFR);

%% baseline correction
figure
cfg = [];
cfg.baseline     = [-2 0];
cfg.baselinetype = 'relative';
cfg.channel      = 'R016-2012-10-03-CSC04a';
ft_singleplotTFR(cfg, TFR);

%% Frequency-dependent windowing
cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.channel      = 'R016-2012-10-03-CSC04a';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 10:2:100; % frequencies of interest: 10hz to 100hz in steps of 2Hz(1/windowSizeinSec)
cfg.t_ftimwin = 20./cfg.foi; 
cfg.toi          = -0.5:0.05:3.5; % times of interest: time window "slides" 
                                  % from -0.5 to 3.5 sec in steps of 0.05 sec (50 ms)
 
TFR = ft_freqanalysis(cfg, data_trl);
 
figure
cfg = []; 
cfg.baseline     = [-2 0];
cfg.baselinetype = 'relative';
cfg.channel = 'R016-2012-10-03-CSC04a';
ft_singleplotTFR(cfg, TFR);

%% statistical tests
% a statistical comparison of the pre-nosepoke baseline and the oscillation 
% patterns following the nosepoke

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'R016-2012-10-03-CSC04a';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100;
cfg.keeptrials   = 'yes'; % need this for stats later
cfg.t_ftimwin    = 20./cfg.foi;  % 20 cycles per time window
 
cfg.toi          = -1:0.05:0; % pre-nosepoke baseline
TFR_pre = ft_freqanalysis(cfg, data_trl);
 
cfg.toi          = 0:0.05:1; % post-nosepoke
TFR_post = ft_freqanalysis(cfg, data_trl);
 
TFR_pre.time = TFR_post.time; % time axes should be identical for comparison
%% t-test
cfg = [];
cfg.channel     = 'R016-2012-10-03-CSC04a';
cfg.latency     = 'all';
cfg.trials      = 'all';
cfg.frequency   = 'all';
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.avgoverfreq = 'no';
cfg.parameter   = 'powspctrm';
cfg.method      = 'stats';
cfg.statistic   = 'ttest2';
cfg.alpha       = 0.05;
 
nTrials1 = size(TFR_pre.powspctrm,1); 
nTrials2 = size(TFR_post.powspctrm,1);
 
cfg.design = cat(2,ones(1,nTrials1),2*ones(1,nTrials2)); % two conditions
cfg.ivar = 1; % dimension of design var which contains the independent variable (group)
 
stat = ft_freqstatistics(cfg,TFR_post,TFR_pre);
 
cfg.parameter = 'stat';
ft_singleplotTFR(cfg,stat); % plot the t-statistic

%% multichannel data
%  load the data
fc = FindFiles('*.ncs'); % get filenames of all LFPs recorded
data_all = ft_read_neuralynx_interp(fc); % load them all -- this will take a while
data_all.hdr.Fs = data_all.fsample; % for some reason this is missing from the header
 
% define layout for later plotting
cfg = [];
cfg.layout = 'ordered'; cfg.channel = data.label;
layout = ft_prepare_layout(cfg, data_all);

%  trial definition 
cfg = [];
cfg.t = cat(2,getd(evt,'n0'),getd(evt,'n1'));
cfg.mode = 'nlx';
cfg.hdr = data_all.hdr;
cfg.twin = [-1 4];
 
trl = ft_maketrl(cfg);
 
cfg = [];
cfg.trl = trl;
data_trl = ft_redefinetrial(cfg,data_all); 
 
% spectrogram computations 
cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:100; % frequencies of interest
cfg.toi          = -0.5:0.05:3.5; % times of interest
cfg.t_ftimwin = 20./cfg.foi;
 
TFR = ft_freqanalysis(cfg, data_trl);

% final plot: baseline-corrected, event-aligned spectrograms 
% for all 16 channels in this session
% There is also an average shown (the subplot on the lower right). 
figure
cfg = [];
cfg.baseline     = [-2 0];
cfg.baselinetype = 'relative';
cfg.layout = layout;
 
ft_multiplotTFR(cfg, TFR); % note this is now multiplot rather than singleplot





