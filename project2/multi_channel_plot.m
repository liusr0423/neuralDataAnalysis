cfg = [];
 
cfg.eventList = {'Feeder 0 nosepoke','Feeder 1 nosepoke', ...
    '1 pellet cue','3 pellet cue','5 pellet cue', ...
    '1 pellet dispensed','3 pellet dispensed','5 pellet dispensed'};
 
cfg.eventLabel = {'n0','n1', ...
    'c1','c3','c5', ...
    'd1','d3','d5'};
 
evt = LoadEvents(cfg);

%% load the data
fc = FindFiles('*.ncs'); % get filenames of all LFPs recorded
data_all = ft_read_neuralynx_interp(fc); % load them all -- this will take a while
data_all.hdr.Fs = data_all.fsample; % for some reason this is missing from the header
 
%% define layout for later plotting
cfg = [];
cfg.layout = 'ordered'; 
cfg.channel = data_all.label;
layout = ft_prepare_layout(cfg, data_all);

%%
cfg = [];
cfg.t = cat(2,getd(evt,'d1'),getd(evt,'d3'),getd(evt,'d5'));
cfg.mode = 'nlx';
cfg.hdr = data_all.hdr;
cfg.twin = [-1 4];
 
trl = ft_maketrl(cfg);
 
cfg = [];
cfg.trl = trl;
data_trl = ft_redefinetrial(cfg,data_all); 
 
%%
cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:100; % frequencies of interest
cfg.toi          = -0.5:0.05:3.5; % times of interest
cfg.t_ftimwin = 20./cfg.foi;
 
TFR = ft_freqanalysis(cfg, data_trl);

figure
cfg = [];
cfg.baseline     = [-2 0];
cfg.baselinetype = 'relative';
cfg.layout = layout;
colormap jet
ft_multiplotTFR(cfg, TFR); % note this is now multiplot rather than singleplot