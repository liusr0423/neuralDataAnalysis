%% some initial path settings
restoredefaultpath; % start with a clean slate
clear classes;

cd('/Users/liusirui/Documents/MATLAB/fieldtrip'); % replace paths with yours
ft_defaults;

cd('/Users/liusirui/Documents/MATLAB/neuraldata-w16/shared'); % replace paths with yours
p = genpath(pwd); % create list of all folders from here
addpath(p);

set(0,'DefaultAxesFontSize',10)

%% define and cd to data folder
df = 'R020';
dataDir = ['/Users/liusirui/Documents/MATLAB/class/neuralDataAnalysis/Data/',df]; % replace data paths with yours
cd(dataDir);
dataFiles = dir([df,'*']);
nSessions = length(dataFiles);

%% for each data session
for ss = 1:nSessions  % loop over sessions
    this_session = dataFiles(ss).name;
    cd(this_session);
    LoadExpKeys; % load current data session 
    ExpKeys
    
    % find LFP file recorded from Hippocampus
    HC_idx = find(ExpKeys.TetrodeTargets == strmatch('Hippocampus',ExpKeys.Target));
    HC_df = {};
    for ii = 1:length(HC_idx)
        HC_df{ii}(:) = dir([this_session(1:end-length('promoted_')),sprintf('-CSC0%d*.ncs',HC_idx(ii))]);
    end
    
    %% Overall comparison of HC-Str coherence
    
    % Load two hippocampus and ventral striatum LFPs
    cfg = [];
    cfg.fc = {HC_df{1}(1).name,ExpKeys.goodGamma{1}};
    cfg.label = {'HC','Str'};
    csc = LoadCSC(cfg);
    csc = restrict(csc,ExpKeys.TimeOnTrack(1),ExpKeys.TimeOffTrack(2)); % restrict to task
    
    % Compute PSD for each LFP and the overall coherence between the pair
    Fs = csc.cfg.hdr{1}.SamplingFrequency;
    wsize = 2048;
    nS = length(csc.label);
    P = {};F = {}; C = {}; Fc = {};
    for iS = 1:nS
        [P{iS},F{iS}] = pwelch(getd(csc,csc.label{iS}),hanning(wsize),wsize/2,2*wsize,Fs);
        for iS2 = iS+1:nS
            [C{iS,iS2},Fc{iS}] = mscohere(getd(csc,csc.label{iS}),getd(csc,csc.label{iS2}),hanning(wsize),wsize/2,2*wsize,Fs);
        end
    end
    
    % plot
    figure
    subplot(121)
    cols = 'kgm';
    for iS = 1:nS
        h(iS) = plot(F{iS},10*log10(P{iS}),cols(iS),'LineWidth',2); hold on;
    end
    set(gca,'XLim',[0 100],'XTick',0:10:100,'FontSize',12); grid on;
    legend(h,csc.label,'Location','Northeast'); legend boxoff;
    xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    title(sprintf('%s',this_session(1:end-9)));
    
    subplot(122); clear h;
    h(1) = plot(Fc{1},C{1,2},'LineWidth',2); hold on;
    set(gca,'XLim',[0 100],'XTick',0:10:100,'FontSize',12); grid on;
    legend(h,{'HC-vStr'},'Location','Northeast'); legend boxoff;
    xlabel('Frequency (Hz)'); ylabel('Coherence');
    title(sprintf('%s',this_session(1:end-9)));
    % the PSD plot on the left shows that hippocampus has a clear theta
    % (verified the correct hippocampus file selection) and ventral striatum
    % has a large gamma component;
    % coherence plot shows an overall HC-vStr coherence peak at 
    % around theta band
    
    %% Task-based coherence analysis using fieldtrip
    % load the data
    cfg = [];
    cfg.fc = cat(2,HC_df{1}(1).name,ExpKeys.goodGamma(1));
    data = ft_read_neuralynx_interp(cfg.fc);
    data.label = {'HC','Str'};
    data.hdr.Fs = data.fsample;
    
    % trilify: trials when rat nosepoked into the reward receptacle
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
    %%
    % spectral decomposition
    cfg              = [];
    cfg.output       = 'powandcsd';  % output power and cross-spectral density
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 1:1:100;      % frequencies to use
    cfg.t_ftimwin    = 20./cfg.foi;  % frequency-dependent, 20 cycles per time window
    cfg.keeptrials   = 'yes';
    cfg.channel      = {'HC' 'Str'};
    cfg.channelcmb   = {'HC' 'Str'}; % channel pairs to compute csd for
    cfg.toi          = -2:0.05:4.5;
    TFR_conn         = ft_freqanalysis(cfg, data_trl);
    %%
    % compute and plot pairwise phase consistency(ppc) as a function of 
    % time and frequency
    cfg            = [];
    cfg.method     = 'ppc';
    fd_ppc         = ft_connectivityanalysis(cfg,TFR_conn);

    figure;
    lbl= cat(2,fd_ppc.labelcmb{1},'-',fd_ppc.labelcmb{2}); % get the label of this pair
    imagesc(fd_ppc.time,fd_ppc.freq,sq(fd_ppc.ppcspctrm(1,:,:))); 
    axis xy; colorbar; colormap jet;
    xlabel('time (s)'); ylabel('Frequency (Hz)'); title(lbl);
    % time-frequency coherence spectrums for all sessions 
    % show similar theta band (~8Hz) coherence peak 1s around
    % nosepoke onset (t=0) and also similar beta band (~15-20Hz) coherence 
    % peak before nosepoke onset 
    %%
    % compute coherence from the cross-spectrum and the indvidual spectra
    cfg            = [];
    cfg.method     = 'coh';
    fd_conn        = ft_connectivityanalysis(cfg,TFR_conn);
    
    % plot coherence spectrum
    figure;
    subplot(211)
    lbl = cat(2,fd_conn.labelcmb{1}, '-', fd_conn.labelcmb{2});
    h = plot(fd_conn.freq, nanmean(sq(fd_conn.cohspctrm(1,:,:)),2));
    hold on;
    xlim([0,100]);set(gca,'XTick',0:10:100);
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    legend(h,lbl);
    title(sprintf('coherence spectrum (%s)',this_session(1:end-9)))
    % coherence spectrums for all sessions show similar coherence peak 
    % around theta band (~8Hz) 
    
    %% Task-based causality analysis: phase-slope index
    % frequency decomposition
    cfg_TFR = [];
    cfg_TFR.output = 'fourier';
    cfg_TFR.method = 'mtmfft';
    cfg_TFR.taper = 'hanning';
    cfg_TFR.foi = 1:1:100;
    cfg_TFR.channel = {'HC','Str'};
    cfg_TFR.channelcmb = {'HC','Str'};
    TFR_psi = ft_freqanalysis(cfg_TFR,data_trl);
    
    % compute psi
    cfg_psi = [];
    cfg_psi.method = 'psi';
    cfg_psi.bandwidth = 4; % number of frequencies to compute slope over
    cfg_psi.channel = {'HC','Str'};
    cfg_psi.channelcmb = {'HC','Str'};   
    psi = ft_connectivityanalysis(cfg_psi,TFR_psi);

    % plot psi
    subplot(212)
    hold on
    plot(psi.freq,sq(psi.psispctrm(1,2,:)).*(360/(2*pi)^2)); % convert to degrees/Hz
    plot(psi.freq, psi.freq*0, '--r')
    title(sprintf('phase-slope (%s)',this_session(1:end-9)));
    xlim([0 100]);
    set(gca,'XTick',0:10:100);
    xlabel('Frequency (Hz)'); ylabel('Phase slope (?/Hz)');
    % positive means HC leads Str
    
    % five out of seven psi plots show similar positive phase slope around
    % theta band (~8Hz) (HC->Str) and the rest is negative around theta band
    % (Str->HC)
    %% store everything for this session
    All{ss}.psi.freq = psi.freq;
    All{ss}.psi.psispctrm = sq(psi.psispctrm(1,2,:)).*(360/(2*pi)^2);
    All{ss}.coh.cohspctrm = nanmean(sq(fd_conn.cohspctrm(1,:,:)),2);
    
    cd(dataDir); % back to root data folder
end
