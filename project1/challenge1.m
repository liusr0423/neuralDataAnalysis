%% load the data files
rootDir = '/Users/liusirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R016-2012-10-03';
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
%% plot the spectrogram (frequency-dependent windowing) for each event

for n = 1:length(evt.label)
    cfg = [];
    cfg.t = getd(evt,evt.label{n});
    cfg.mode = 'nlx';
    cfg.hdr = data.hdr;
    cfg.twin = [-1 4];
    
    trl = ft_maketrl(cfg);
    
    cfg = [];
    cfg.trl = trl;
    data_trl = ft_redefinetrial(cfg,data);
    data_trl = ft_removeMissingDataTrials([],data_trl);

    cfg              = []; % start with empty cfg
    cfg.output       = 'pow';
    cfg.channel      = 'R016-2012-10-03-CSC04a';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 10:2:100; % frequencies of interest: 10hz to 100hz in steps of 1Hz(1/windowSizeinSec)
    cfg.t_ftimwin    = 20./cfg.foi;
    cfg.toi          = -.5:0.05:3.5; % times of interest: time window "slides"
                                      % from -0.5 to 3.5 sec in steps of 0.05 sec (50 ms)
    
    TFR = ft_freqanalysis(cfg, data_trl);
    
    figure(n)
    cfg = [];
    cfg.baseline     = [-2 0];
    cfg.baselinetype = 'relative';
    cfg.channel = 'R016-2012-10-03-CSC04a';
    ft_singleplotTFR(cfg, TFR);
    title(evt.label{n},'FontSize',20);
    colormap jet
end

%% plot event-triggered LFP traces for each event
cd('/Users/liusirui/Documents/MATLAB/class/neuralDataAnalysis/project1');
cfg = [];
cfg.eventLabel = evt.label;
cfg.eventTimes = evt.t(5);
cfg.twin = [-.5,3.5];
cfg.f = [70,90];
eventLFPplot(cfg,csc)







