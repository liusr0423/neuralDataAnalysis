%% load the data files
rootDir = '/Users/sirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R016-2012-10-03';
cd(rootDir);
cfg = [];
cfg.fc = {'R016-2012-10-03-CSC04a.ncs'};
csc = LoadCSC(cfg);

%% Load trials info
cfg = [];
cfg.eventList = {'Feeder 0 nosepoke','Feeder 1 nosepoke', ...
    '1 pellet cue','3 pellet cue','5 pellet cue', ...
    '1 pellet dispensed','3 pellet dispensed','5 pellet dispensed'};
cfg.eventLabel = {'n0','n1', ...
    'c1','c3','c5', ...
    'd1','d3','d5'}; 

evt = LoadEvents(cfg);

%% Loading Neuralynx data into FieldTrip
fc = {'R016-2012-10-03-CSC04a.ncs'};
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
data_trl  = ft_removeMissingDataTrials([],data_trl);
%% Frequency-dependent windowing
cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.channel      = 'R016-2012-10-03-CSC04a';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 10:2:100; % frequencies of interest: 10hz to 100hz in steps of 2Hz(1/windowSizeinSec)
cfg.t_ftimwin    = 20./cfg.foi; 
cfg.toi          = -0.5:0.05:3.5; % times of interest: time window "slides" 
                                  % from -0.5 to 3.5 sec in steps of 0.05 sec (50 ms)
 
TFR = ft_freqanalysis(cfg, data_trl);
 
figure
cfg = []; 
cfg.baseline     = [-2 0];
cfg.baselinetype = 'relative';
cfg.channel = 'R016-2012-10-03-CSC04a';
ft_singleplotTFR(cfg, TFR);
%% plot event-triggered LFP traces for each event
cfg = [];
cfg.eventTimes = evt.t;
cfg.twin = [-1,3];
cfg.f = [50,60];
eventLFPplot(cfg,csc)







