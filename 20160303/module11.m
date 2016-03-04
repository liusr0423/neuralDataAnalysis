%% Diversion: the Wiener-Khinchin theorem
% the fourier transform of a given signal's autocorrelation corresponds to
% the fourier transform of the signal itself
Fs = 500; dt = 1./Fs;
t = [0 2]; tvec = t(1):dt:t(2)-dt;
 
f1 = 8;
data1 = sin(2*pi*f1*tvec)+0.1*randn(size(tvec));
 
[acf,lags] = xcorr(data1,100,'coeff'); 
% xcorr(A),is the auto-correlation sequence
% 'coeff' = normalizes the sequence so that the auto-correlations at zero
% lag are identically 1.0
lags = lags.*(1./Fs); % convert samples to time
plot(lags,acf); grid on;

%% cross-correlation between two signals
f2 = 8;
data2 = sin(2*pi*f2*tvec+pi/4)+0.1*randn(size(tvec)); % phase-shifted version of data1
 
[ccf,lags] = xcorr(data1,data2,100,'coeff'); % now a cross-correlation
lags = lags.*(1./Fs); % convert samples to time
plot(lags,ccf); grid on;
xlabel('lags(s)');ylabel('xcorr');title('data1-data2');
% the plot shows the probablity of signal 2 at various time lags (between 
% -.2 and .2s) given that signal 1 at time 0

%% coherence
figure;
subplot(221);
plot(tvec,data1,'r',tvec,data2,'b'); legend({'signal 1','signal 2'});
title('raw signals');
 
[Pxx,F] = pwelch(data1,hanning(250),125,length(data1),Fs);
[Pyy,F] = pwelch(data2,hanning(250),125,length(data1),Fs);
subplot(222)
plot(F,abs(Pxx),'r',F,abs(Pyy),'b'); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('PSD');
 
[Pxy,F] = cpsd(data1,data2,hanning(250),125,length(data1),Fs);
% cross-spectrum function cpsd(data1,data2,window,overlaap,nFFT,Fs)
subplot(223)
plot(F,abs(Pxy)); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('cross-spectrum');
 
[acf,lags] = xcorr(data1,data2,100,'coeff');
lags = lags.*(1./Fs); % convert samples to time
 
subplot(224)
plot(lags,acf); grid on;
xlabel('time lag (s)'); ylabel('correlation ({\itr})'); title('xcorr');
%% if data1 twice as large the amplitude
data1 = 2 * (sin(2*pi*f1*tvec)+0.1*randn(size(tvec)));
figure;
subplot(221);
plot(tvec,data1,'r',tvec,data2,'b'); legend({'signal 1','signal 2'});
title('raw signals');
 
[Pxx,F] = pwelch(data1,hanning(250),125,length(data1),Fs);
[Pyy,F] = pwelch(data2,hanning(250),125,length(data1),Fs);
subplot(222)
plot(F,abs(Pxx),'r',F,abs(Pyy),'b'); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('PSD');
 
[Pxy,F] = cpsd(data1,data2,hanning(250),125,length(data1),Fs);
% cross-spectrum function cpsd(data1,data2,window,overlaap,nFFT,Fs)
subplot(223)
plot(F,abs(Pxy)); xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('power'); title('cross-spectrum');
 
[acf,lags] = xcorr(data1,data2,100,'coeff');
lags = lags.*(1./Fs); % convert samples to time
 
subplot(224)
plot(lags,acf); grid on;
xlabel('time lag (s)'); ylabel('correlation ({\itr})'); title('xcorr');
%  the cross-spectrum depends on the amplitude of the input signals.
% Thus, coherence normalizes the cross-spectrum by the spectra of the 
% individual signals.

%% compute the coherence by normalizing the cross-spectrum
% Verify that now there is no change in coherence when scaling data1 by a 
% factor 2 as above.
coeff = 2;
data1 = sin(2*pi*f1*tvec)+0.1*randn(size(tvec));
data1_scaled = coeff * data1;
data2 = sin(2*pi*f2*tvec+pi/4)+0.1*randn(size(tvec)); % phase-shifted data1
data1 = {data1,data1_scaled};
C =cell(1,2);
for n = 1:2
[Pxx,F] = pwelch(data1{n},hanning(250),125,length(data1{n}),Fs);
[Pyy,F] = pwelch(data2,hanning(250),125,length(data1{n}),Fs);
[Pxy,F] = cpsd(data1{n},data2,hanning(250),125,length(data1{n}),Fs);
C{n} = (abs(Pxy).^2)./(Pxx.*Pyy);
end

tf = isequal(C{1},C{2})
figure

plot(F,C{1},'b');
hold on
plot(F,C{2},'r');
set(gca,'Xlim',[0 20],'XTick',0:2:20)
title('Magnitude-squared Coherence')
xlablel('Frequency (Hz)')

%% cross spectrum phase
data1 = sin(2*pi*f1*tvec)+0.1*randn(size(tvec));
data2 = sin(2*pi*f2*tvec+pi/4)+0.1*randn(size(tvec)); % phase-shifted data1
[Pxy,F] = cpsd(data1,data2,hanning(250),125,length(data1),Fs);
figure
plot(F,angle(Pxy));
xlabel('Frequency (Hz)')
ylabel('Lag (\times\pi rad)')
title('Cross Spectrum Phase')
% cross-spectrum estimates are spaced at 500/1000(Fs/nFFT)=.5 so to return
% the phase estimates at 8Hz which should be close to the true values
phi8 = -angle(Pxy(8/.5+1)); % first freq bin corresponds to 0
lag8 = phi8/pi;% phi8 should be close to +pi/4

%% properties of the coherence measure: just verify some cases where we break the phase relationship
% e.g.1.two signals that have an 8Hz component in their power spectra are not necessarily coherent:
f = 0.5; % freq modulation (Hz) 
f2 = 8;f1 = 8;
m = 4; % freq modulation strength
wsz = 250; % window size 
t = [0 2]; tvec = t(1):dt:t(2)-dt; %  The coherence estimate should clean up somewhat as you increase the data length
data1 = sin(2*pi*f1*tvec)+0.1*randn(size(tvec));
data2 = sin(2*pi*f2*tvec+pi/4)+0.1*randn(size(tvec)); 

figure
subplot(421)
s2 = data2;
plot(tvec,s2,tvec,data1); title('signal 1 - constant phase');
 
subplot(422)
s3 = sin(2*pi*f2*tvec + m.*sin(2*pi*f*tvec - pi/2)) + 0.1*randn(size(tvec));
plot(tvec,s3,tvec,data1); title('signal 2 - varying phase');
 
subplot(423)
[Ps2,F] = pwelch(s2,hanning(wsz),wsz/2,length(data2),Fs);
plot(F,abs(Ps2)); title('PSD');
 
subplot(424)
[Ps3,F] = pwelch(s3,hanning(wsz),wsz/2,length(data2),Fs);
plot(F,abs(Ps3)); title('PSD');
 
subplot(425)
[C,F] = mscohere(data1,s2,hanning(wsz),wsz/2,length(data1),Fs); % shortcut to obtain coherence
plot(F,C); title('coherence'); xlabel('Frequency (Hz)');
 
subplot(426)
[C,F] = mscohere(data1,s3,hanning(wsz),wsz/2,length(data1),Fs);
plot(F,C); title('coherence'); xlabel('Frequency (Hz)');
 
[acf,lags] = xcorr(data1,s2,100,'coeff');
lags = lags.*(1./Fs); % convert samples to time
 
subplot(427)
plot(lags,acf); grid on;
xlabel('time lag (s)'); ylabel('correlation ({\itr})'); title('xcorr');
 
[acf,lags] = xcorr(data1,s3,100,'coeff');
lags = lags.*(1./Fs); % convert samples to time
 
subplot(428)
plot(lags,acf); grid on;
xlabel('time lag (s)'); ylabel('correlation ({\itr})'); title('xcorr');
%% e.g.2: 
% two signals both have 40hz components in PSD,but the times at which the 
% 40Hz oscillation is present do not actually overlap between the two signals
wsize = 500;
 
Fs = 500; dt = 1./Fs;
t = [0 2];
 
tvec = t(1):dt:t(2)-dt;
f1 = 40; f2 = 40;
 
% generate some square waves
mod1 = square(2*pi*4*tvec,20); mod1(mod1 < 0) = 0;
mod2 = square(2*pi*4*tvec+pi,20); mod2(mod2 < 0) = 0;
 
data1 = sin(2*pi*f1*tvec); data1 = data1.*mod1 + 0.01*randn(size(tvec));
data2 = sin(2*pi*f2*tvec); data2 = data2.*mod2 + 0.01*randn(size(tvec)) ;
figure
subplot(221);
plot(tvec,data1,'r',tvec,data2,'b'); legend({'signal 1','signal 2'});
title('raw signals');
 
[P1,F] = pwelch(data1,hanning(wsize),wsize/2,length(data2),Fs);
[P2,F] = pwelch(data2,hanning(wsize),wsize/2,length(data2),Fs);
subplot(222)
plot(F,abs(P1),'r',F,abs(P2),'b'); title('PSD');
 
subplot(223);
[C,F] = mscohere(data1,data2,hanning(wsize),wsize/2,length(data1),Fs);
plot(F,C); title('coherence'); xlabel('Frequency (Hz)');
 
[ccf,lags] = xcorr(data1,data2,100,'coeff');
lags = lags.*(1./Fs); % convert samples to time
 
subplot(224)
plot(lags,ccf); grid on;
xlabel('time lag (s)'); ylabel('correlation ({\itr})'); title('xcorr');

%% Application to real data
% load three simultaneously recorded LFPs, two from the same structure (but 
% a different electrode, both in ventral striatum) and one from a different 
% but anatomically related structure (hippocampus)
clear all
cd ('/Users/sirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R016-2012-10-03')
LoadExpKeys;
cfg = []; 
cfg.fc = cat(2,ExpKeys.goodGamma(1:2),ExpKeys.goodTheta(1)); % cell array containing filenames to load
cfg.label = {'vStr1','vStr2','HC'};
csc = LoadCSC(cfg);
 
csc = restrict(csc,ExpKeys.TimeOnTrack(1),ExpKeys.TimeOffTrack(2)); % restrict to task
%% compute PSDs for each signal & coherence between signal pairs of interest

Fs = csc.cfg.hdr{1}.SamplingFrequency; 
wsize = 2048;
nS = length(csc.label);
for iS = 1:nS
    [P{iS},F{iS}] = pwelch(getd(csc,csc.label{iS}),hanning(wsize),wsize/2,2*wsize,Fs);
 
    for iS2 = iS+1:nS
        [C{iS,iS2},Fc{iS}] = mscohere(getd(csc,csc.label{iS}),getd(csc,csc.label{iS2}),hanning(wsize),wsize/2,2*wsize,Fs);
    end
end

%% plot
subplot(121)
cols = 'kgm';
for iS = 1:nS
h(iS) = plot(F{iS},10*log10(P{iS}),cols(iS),'LineWidth',2); hold on;
end
set(gca,'XLim',[0 150],'XTick',0:25:150,'FontSize',12); grid on;
legend(h,csc.label,'Location','Northeast'); legend boxoff;
xlabel('Frequency (Hz)'); ylabel('Power (dB)'); 
 
subplot(122); clear h;
h(1) = plot(Fc{1},C{1,2},'LineWidth',2); hold on;
h(2) = plot(Fc{1},C{1,3},'r','LineWidth',2);
set(gca,'XLim',[0 150],'XTick',0:25:150,'FontSize',12); grid on;
legend(h,{'vStr1-vStr2','vStr1-HC'},'Location','Northeast'); legend boxoff;
xlabel('Frequency (Hz)'); ylabel('Coherence');



%% Comparison of vStr-HC coherence between experimental conditions
% is the coherence between hippocampus and ventral striatum modulated by
% task events - if there is a change in coherence between approach to the 
% reward site and reward receipt. 
% using fieldtrip
cd('/Users/sirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R016-2012-10-03'); 
LoadExpKeys;
 
cfg.fc = cat(2,ExpKeys.goodGamma(1:2),ExpKeys.goodTheta(1));
data = ft_read_neuralynx_interp(cfg.fc);
data.label = {'vStr1','vStr2','HC'};


%% trialify: extract task-specfic timestamps
% for this dataset toi are the times rat nosepoke into the reward
% receptacles in anticipation of receiving a number of pellets

data.hdr.Fs = data.fsample;
 
cfg = [];
cfg.trialfun = 'ft_trialfun_lineartracktone2';
cfg.trialdef.hdr = data.hdr;
cfg.trialdef.pre = 2.5; 
cfg.trialdef.post = 5; % define time window of interest
 
cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'; this and what follows are all task-specific
cfg.trialdef.location = 'both'; % could be 'left', 'right', 'both'
cfg.trialdef.block = 'both'; % could be 'value', 'risk', 'both'
cfg.trialdef.cue = {'c1','c3','c5'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
 
[trl, event] = ft_trialfun_lineartracktone2(cfg);
cfg.trl = trl;
 
data_trl = ft_redefinetrial(cfg,data);

%% compute the trial-averaged cross-spectrum
% similar to the code used for computing spectrograms but changed
% cfg.output from 'pow' to 'powandcsd' (csd is used for cross-spectral
% density)
cfg              = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100; % frequencies to use
cfg.t_ftimwin    = 20./cfg.foi;  % frequency-dependent, 20 cycles per time window
cfg.keeptrials   = 'yes';
cfg.channel      = {'vStr1','vStr2', 'HC'};
cfg.channelcmb   = {'vStr2', 'HC'; 'vStr2', 'vStr1'};  % channel pairs to compute csd for

cfg.toi          = -2:0.05:0; % toi: pre-nosepoke baseline (time 0 is time of nosepoke)

TFR_pre = ft_freqanalysis(cfg, data_trl);

%% compute the coherence from cross-spectrum and the individual spectra
cfg            = [];
cfg.method     = 'coh'; % compute coherence; other measures of connectivity are also available
fd             = ft_connectivityanalysis(cfg,TFR_pre);
%% plot
figure(1);
cols = 'rgb';
for iCmb = 1:size(fd.labelcmb,1)
    lbl{iCmb} = cat(2,fd.labelcmb{iCmb,1},'-',fd.labelcmb{iCmb,2});
    temp = nanmean(sq(fd.cohspctrm(iCmb,:,:)),2);
    h1(iCmb) = plot(fd.freq,temp,cols(iCmb));
    hold on;
end
xlabel('Frequency(Hz)')
ylabel('Coherence')
hold on;

%% post-nosepoke period
cfg              = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100; % frequencies to use
cfg.t_ftimwin    = 20./cfg.foi;  % frequency-dependent, 20 cycles per time window
cfg.keeptrials   = 'yes';
cfg.channel      = {'vStr1','vStr2', 'HC'};
cfg.channelcmb   = {'vStr2', 'HC'; 'vStr2', 'vStr1'};  % channel pairs to compute csd for

cfg.toi          = 0:0.05:2; % toi: pre-nosepoke baselfine (time 0 is time of nosepoke)

TFR_post = ft_freqanalysis(cfg, data_trl);

cfg            = [];
cfg.method     = 'coh'; % compute coherence; other measures of connectivity are also available
fd             = ft_connectivityanalysis(cfg,TFR_post);

cols = 'bmk';
for iCmb = 1:size(fd.labelcmb,1)
    lbl{iCmb} = cat(2,fd.labelcmb{iCmb,1},'-',fd.labelcmb{iCmb,2});
    temp = nanmean(sq(fd.cohspctrm(iCmb,:,:)),2);
    h2(iCmb) = plot(fd.freq,temp,cols(iCmb));
    hold on;
end
legend([h1,h2],[lbl{1},' pre-nosepoke'],[lbl{2},' pre-nosepoke',lbl{1},' post-nosepoke'],[lbl{2},' post-nosepoke']);

% pre-nosepoke: coherence peak around 15Hz 
% post-nosepoke: coherence peak around 10Hz


%% time-frequency coherence analysis
% coherogram: coherence as a function of time and frequency
iC = 1; % which signal pair to plot
lbl = [fd.labelcmb{1,:}]; % get the label of this pair
imagesc(fd.time,fd.freq,sq(fd.cohspctrm(iC,:,:))); 
axis xy; colorbar; xlabel('time (s)'); ylabel('Frequency (Hz)'); 
title(lbl);colormap jet
% We can improve the robustness of our estimate by giving up some time resolution
% or averaging over more trials
%% ft()connectivityanalysis()
cfg            = [];
cfg.method     = 'ppc';
fd             = ft_connectivityanalysis(cfg,TFR_post);
iC = 1; % which signal pair to plot
lbl = [fd.labelcmb{1,:}]; % get the label of this pair
imagesc(fd.time,fd.freq,sq(fd.ppcspctrm(iC,:,:))); 
axis xy; colorbar; xlabel('time (s)'); ylabel('Frequency (Hz)'); 
title(lbl);colormap jet

%% Granger causality: test directionality
% generate some artificial data: 1000 trials of 5 seconds each of
% independent white noise for two signals X and Y
cfg             = [];
cfg.ntrials     = 1000;
cfg.triallength = 5; % in seconds
cfg.fsample     = 1000;
cfg.nsignal     = 2; % two signals, X and Y, which start out as identical white noise
 
cfg.method      = 'linear_mix';
cfg.mix         = [0; 0]; % multiply white noise for X and Y by this
cfg.delay       = [0; 0]; % Y is n samples delayed relative to X (both 0)
cfg.bpfilter    = 'no';
cfg.absnoise    = 1; % add independent noise to both signals, so now X and Y should be independent
 
data            = ft_connectivitysimulation(cfg);
data.label      = {'X','Y'};

%% Verify that indeed the two signals X and Y are uncorrelated
a = [];
for n = 1:length(data.trial)
temp = corrcoef(data.trial{n}(1,:),data.trial{n}(2,:));a(n)=temp(1,2);
end
hist(a);
%% fit AR model
cfg_ar         = [];
cfg_ar.order   = 3; % how far back to estimate coefficients for
cfg_ar.toolbox = 'bsmart';
mdata          = ft_mvaranalysis(cfg_ar, data);

%% plot the coefficients
% the extent that we can predict each signal separately based on its own 
% past, and then how much that prediction can be improved by knowledge of 
% the other signal.

figure; subplot(221)
 
labels = {'X->X','X->Y';'Y->X','Y->Y'}; cols = 'rgbc';
nP = 0;
for iI = 1:cfg.nsignal   
    for iJ = 1:cfg.nsignal
        nP = nP + 1;
        h(nP) = plot(1:cfg_ar.order,sq(mdata.coeffs(iI,iJ,:)),cols(nP));
        hold on;
        plot(1:cfg_ar.order,sq(mdata.coeffs(iI,iJ,:)),'.','MarkerSize',20,'Color',cols(nP));
 
    end
end
set(gca,'FontSize',18,'LineWidth',1); box off;
set(h,'LineWidth',2);
xlabel('lag (samples)'); ylabel('coefficient');
title('cfg.delay = [0; 0];');
legend(h,labels(:));
% small coefficient values -> uncorrelated signals 
% we cannot predict anything about our signal based on its past ? 
% the definition of white noise

%% X -> Y case 
cfg.mix         = [0.8; 0.8]; % X and Y are identical white noise with amplitude 0.8
cfg.absnoise    = 0.2; % add amplitude 0.2 *independent noise
cfg.delay       = [0; 2]; % advance Y 2 samples relative to X
 
data            = ft_connectivitysimulation(cfg);
data.label      = {'X','Y'};

cfg_ar         = [];
cfg_ar.order   = 3; % how far back to estimate coefficients for
cfg_ar.toolbox = 'bsmart';
mdata          = ft_mvaranalysis(cfg_ar, data);

labels = {'X->X','X->Y';'Y->X','Y->Y'}; cols = 'rgbc';
nP = 0;
subplot(222)
for iI = 1:cfg.nsignal   
    for iJ = 1:cfg.nsignal
        nP = nP + 1;
        h(nP) = plot(1:cfg_ar.order,sq(mdata.coeffs(iI,iJ,:)),cols(nP));
        hold on;
        plot(1:cfg_ar.order,sq(mdata.coeffs(iI,iJ,:)),'.','MarkerSize',20,'Color',cols(nP));
 
    end
end
set(gca,'FontSize',18,'LineWidth',1); box off;
set(h,'LineWidth',2);
xlabel('lag (samples)'); ylabel('coefficient');
title('cfg.delay = [0; 0];');
legend(h,labels(:));
% Note how for the delay case, we correctly estimate that X can be 
% predicted from Y, at the expected delay of 2 samples.


%% Spectrally resolved Granger causality: 
% measures how much of the power in X, not accounted by X itself, can be
% attributed to Y (fit VAR models in the frequency domain)

% generate some artifical data: X and Y are 50% identical signal with 
% frequency content between 50 and 100 Hz, and 50% independent noise
nTrials = 1000;
 
cfg             = [];
cfg.ntrials     = nTrials;
cfg.triallength = 5;
cfg.fsample     = 1000;
cfg.nsignal     = 2;
 
cfg.method      = 'linear_mix';
cfg.mix         = [0.5; 0.5]; % X and Y are identical white noise with amplitude 0.5
cfg.delay       = [0; 4];% advance Y 4 samples relative to X
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [50 100]; % white noise gets filtered in this frequency band
cfg.absnoise    = 0.5; % add amplitude 0.5 * independent noise to both signals
 
data            = ft_connectivitysimulation(cfg);
data.label      = {'X','Y'};




%% frequency decomposition
cfg_TFR = [];
cfg_TFR.channel = {'X','Y'};
cfg_TFR.channelcmb = {'X' 'Y'};
cfg_TFR.method = 'mtmfft';
cfg_TFR.output = 'fourier';
cfg_TFR.foi = 1:1:150;
cfg_TFR.taper = 'hanning';
 
TFR = ft_freqanalysis(cfg_TFR,data);

%% Granger spectra:
cfg_G = [];
cfg_G.method = 'granger';
cfg_G.channel = {'X','Y'};
cfg_G.channelcmb = {'X' 'Y'};
 
C = ft_connectivityanalysis(cfg_G,TFR);
%% plot
%  how much power in X (or in Y) can be predicted based on itself or the
% other signal
figure;
for iP = 1:4
    subplot(2,2,iP);
    plot(C.freq,C.grangerspctrm(iP,:));
    set(gca,'FontSize',14,'YLim',[0 0.5]);
    title([C.labelcmb{iP,:}]);
end
% C.grangerspctrm: x->x,y->x,x->y,y->y
% top right plot (Y->X) has higher coefficients than the reverse (X->Y),
% consistnet with the 4-sample advancement of Y relative to X




%% new data: 
% X is the same signal as Y but twice as large, there is no delay them 
cfg             = [];
cfg.ntrials     = nTrials;
cfg.triallength = 5;
cfg.fsample     = 1000;
cfg.nsignal     = 2;
 
cfg.method      = 'linear_mix';
cfg.mix         = [1; 0.5]; % X bigger than Y
cfg.delay       = [0; 0];
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [50 100]; % white noise gets filtered in this frequency band
cfg.absnoise    = 0.5; % add independent noise to both signals
 
data            = ft_connectivitysimulation(cfg);
data.label      = {'X','Y'};
    

% redo the Granger spectra: X Granger-causes Y 
% this is the case when two (near-identical) signals have different
% signal-to-noise ratio, Granger causality can be easily fooled by this
% solution: to reverse both signals and test again, if Granger asymmetry
% persists after this, we have a tell-tale of a signal-to-noise Granger
% artifact
%% reverse signal
for n = 1:length(a)
    data.trial{n} = flip(data.trial{n},2);
end
% true Granger-causation: spectra should be reversed too




%% Phase-slope index
% phase difference (x-y) as a function of frequency
% if positive slope, x leads y 
% raw phase diffrences (estimating the phase(angle) of the cross-spectrum)
% -> derivative= raw phase slope -> phase slope index 
% (normalized by dividing the raw phase slope at each frequency by its
% standard deviation (estimated using a bootstrap)

% example: generate two signals (Y leads X BUT different signal noise ratio)

nTrials = 1000;
 
cfg             = [];
cfg.ntrials     = nTrials;
cfg.triallength = 5;
cfg.fsample     = 1000;
cfg.nsignal     = 2;
 
cfg.method      = 'linear_mix';
cfg.mix         = [1; 0.3]; % X bigger than Y
cfg.delay       = [0; 4]; % advance Y samples relative to X
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [50 100]; % white noise gets filtered in low gamma band
cfg.absnoise    = 0.5; % add independent noise to both signals
 
data            = ft_connectivitysimulation(cfg);
data.label      = {'X','Y'};

%% compute phase-slope:
% 1. Fourier decomposition
cfg_TFR = [];
cfg_TFR.channel = {'X','Y'};
cfg_TFR.channelcmb = {'X' 'Y'};
cfg_TFR.method = 'mtmfft';
cfg_TFR.output = 'fourier';
cfg_TFR.foi = 1:1:150;
cfg_TFR.taper = 'hanning';
 
TFR = ft_freqanalysis(cfg_TFR,data);

%% 2. use 'psi' as method for connectivity analysis
cfg_psi = [];
cfg_psi.method = 'psi';
cfg_psi.bandwidth = 8; % number of frequencies to compute slope over
cfg_psi.channel = {'X','Y'};
cfg_psi.channelcmb = {'X' 'Y'};
 
C = ft_connectivityanalysis(cfg_psi,TFR);


%% plot the result
figure;
plot(C.freq,sq(C.psispctrm(2,1,:)));
xlabel('Frequency'); ylabel('Phase slope');

% C.psispctrm (row1: x->x,x->y;row2:y->x,y->y)




