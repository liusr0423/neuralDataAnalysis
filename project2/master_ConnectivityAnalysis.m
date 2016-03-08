%% some initial path settings
restoredefaultpath; % start with a clean slate
clear classes; 

cd('/Users/liusirui/Documents/MATLAB/fieldtrip'); % replace paths with yours
ft_defaults;

cd('/Users/liusirui/Documents/MATLAB/neuraldata-w16/shared'); % replace paths with yours
p = genpath(pwd); % create list of all folders from here
addpath(p);

set(0,'DefaultAxesFontSize',18)

%% cd to data folder
df = 'R020';
dataDir = ['../../class/neuralDataAnalysis/Data/',df]; % replace data paths with yours
cd(dataDir);
%% for each data session
dataFiles = dir([df,'*']);
numSessions = length(dataFiles);
this_session = dataFiles(3).name;
cd(this_session);
LoadExpKeys;
ExpKeys
%% find LFP file recorded from Hippocampus 
HC_idx = find(ExpKeys.TetrodeTargets == strmatch('Hippocampus',ExpKeys.Target));
HC_df = {};
for ii = 1:length(HC_idx)
   HC_df{ii}(:) = dir([this_session(1:end-length('promoted_')),sprintf('-CSC0%d*.ncs',HC_idx(ii))]);
end

%% Load two hippocampus and ventral striatum LFPs
cfg = [];
cfg.fc = {HC_df{1}(1).name,ExpKeys.goodGamma{1}};
cfg.label = {'HC','Str'};
csc = LoadCSC(cfg);
 
csc = restrict(csc,ExpKeys.TimeOnTrack(1),ExpKeys.TimeOffTrack(2)); % restrict to task
%% Compute the individual PSD for each LFP and the overall coherence between the pair
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
set(gca,'XLim',[0 150],'XTick',0:10:150,'FontSize',12); grid on;
legend(h,csc.label,'Location','Northeast'); legend boxoff;
xlabel('Frequency (Hz)'); ylabel('Power (dB)'); 

subplot(122); clear h;
h(1) = plot(Fc{1},C{1,2},'LineWidth',2); hold on;
set(gca,'XLim',[0 100],'XTick',0:10:100,'FontSize',12); grid on;
legend(h,{'HC-vStr'},'Location','Northeast'); legend boxoff;
xlabel('Frequency (Hz)'); ylabel('Coherence');
% the PSD plot on the left shows that hippocampus has a clear theta  
% (verified the correct hippocampus file selection) and ventral striatum 
% has a large gamma component and;
% there is an overal HC-vStr coherence peak at around theta band 

%% task-based coherence analysis using ft
cfg = [];
cfg.fc = cat(2,HC_df{1}(1).name,ExpKeys.goodGamma(1));
data = ft_read_neuralynx_interp(cfg.fc);
data.label = {'HC','Str'};
data.hdr.Fs = data.fsample;

%% trilify: trials when rat noespoked into the reward receptacle
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
%% spectral analysis
cfg              = [];
cfg.output       = 'fourier';
%cfg.output       = 'powandcsd';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100; % frequencies to use
cfg.keeptrials   = 'yes';
cfg.channel      = {'HC' 'Str'};
cfg.channelcmb   = {'HC' 'Str'}; % channel pairs to compute csd for

TFR = ft_freqanalysis(cfg, data_trl);

%% compute coherence 
cfg            = [];
cfg.method     = 'coh'; 
fd             = ft_connectivityanalysis(cfg,TFR);

%% plot coherence spectra 
figure;
subplot(211)
lbl= cat(2,fd.label{1},'-',fd.label{2});
h = plot(fd.freq,sq(fd.cohspctrm(1,2,:)));
hold on;

set(gca,'XTick',0:10:100);
xlim([1,100]);
xlabel('Frequency (Hz)'); 
ylabel('Coherence');
legend(h,lbl);
title(sprintf('coherence spectra (%s)',this_session(1:end-9)))
%% causality analysis: phase-slope analysis
% % fourier decomposition
% cfg_TFR = [];
% cfg_TFR.channel      = {'HC','Str'};
% cfg_TFR.channelcmb   = {'HC','Str'}; % channel pairs to compute csd for
% cfg_TFR.method = 'mtmfft';
% cfg_TFR.output = 'fourier';
% cfg_TFR.foi = 1:1:100;
% cfg_TFR.taper = 'hanning';
% TFR = ft_freqanalysis(cfg_TFR, data_trl);

%% compute psi
cfg_psi = [];
cfg_psi.method = 'psi';
cfg_psi.bandwidth = 4; % number of frequencies to compute slope over
cfg_psi.channel = {'HC','Str'};
cfg_psi.channelcmb = {'HC','Str'};
 
C = ft_connectivityanalysis(cfg_psi,TFR);


%% plot psi
subplot(212)
hold on
plot(C.freq,sq(C.psispctrm(1,2,:)));
plot(C.freq, C.freq*0, '--r')
title(sprintf('phase-slope (%s)',this_session(1:end-9)));
set(gca,'XTick',0:10:100);
xlim([0 100]);
xlabel('Frequency (Hz)'); ylabel('Phase slope');
%%
All{1}.psi.freq = C.freq;
All{1}.psi.psispctrm = sq(C.psispctrm(1,2,:)).*(360/(2*pi));
All{1}.coh.cohspctrm = sq(fd.cohspctrm(1,2,:));
% positive means HC leads Str
