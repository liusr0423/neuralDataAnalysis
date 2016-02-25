%% some initial path settings
restoredefaultpath; % start with a clean slate
clear classes; 

cd('/Users/sirui/Documents/MATLAB/fieldtrip'); % replace paths with yours
ft_defaults;

cd('/Users/sirui/Documents/MATLAB/neuraldata-w16/shared'); % replace paths with yours
p = genpath(pwd); % create list of all folders from here
addpath(p);

set(0,'DefaultAxesFontSize',18)

%% load some data files from session R016-2012-10-03-CSC04a
dataDir = '../../class/neuralDataAnalysis/Data/R016-2012-10-03'; % replace data paths with yours
cd(dataDir);
cfg = [];
cfg.fc = {'R016-2012-10-03-CSC04a.ncs'};
csc = LoadCSC(cfg);
Fs = csc.cfg.hdr{1,1}.SamplingFrequency;
%% Load events info for this session
cfg = [];
cfg.eventList = {'Feeder 0 nosepoke','Feeder 1 nosepoke', ...
    '1 pellet cue','3 pellet cue','5 pellet cue', ...
    '1 pellet dispensed','3 pellet dispensed','5 pellet dispensed'};
cfg.eventLabel = {'n0','n1', 'c1','c3','c5','d1','d3','d5'}; 

evt = LoadEvents(cfg);

%% plot the event-triggered spectrogram (frequency-dependent windowing) 
%  for each event: average spectrogram over all trials with 4-second window
% (this may take a while!!)

% Loading Neuralynx data into FieldTrip
fc = {'R016-2012-10-03-CSC04a.ncs'};
data = ft_read_neuralynx_interp(fc);

for n = 1:length(evt.label) % for each trial type
   
    % define and cut the data into trials
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
    
    % construct the event-triggered spectrogram
    cfg              = []; 
    cfg.output       = 'pow';
    cfg.channel      = 'R016-2012-10-03-CSC04a';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 10:2:100; % frequencies of interest: 10hz to 100hz in steps of 2Hz
    cfg.t_ftimwin    = 20./cfg.foi;
    cfg.toi          = -.5:0.05:3.5; % times of interest: time window "slides"
                                     % from -0.5 to 3.5 sec in steps of 0.05 sec (50 ms)
    TFR = ft_freqanalysis(cfg, data_trl);
    
    % plot the spectrograms for each event
    figure(n)
    cfg = [];
    cfg.baseline     = [-2 0];
    cfg.baselinetype = 'relative';
    cfg.channel = 'R016-2012-10-03-CSC04a';
    ft_singleplotTFR(cfg, TFR);
    title(evt.label{n},'FontSize',20);
    set(gca,'FontSize',20);
    set(gcf,'color','w')
    set(gca,'XTick',[-.5:.5:3.5]);
    set(gca,'YTick',[10:10:100]);
    axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');
    colormap jet; colorbar('hide') 
end

%% plot event-triggered trial LFP traces for each event

% cd to project folder
analysisDir = '../../project1'; % replace project path with yours
cd(analysisDir);
%%
% event timestamps and time window
cfg = [];
cfg.eventLabel = {'n0n1','c1c3c5','d1d3d5'};
cfg.eventTimes = {cat(2,getd(evt,'n0'),getd(evt,'n1')),...
    cat(2,getd(evt,'c1'),getd(evt,'c3'),getd(evt,'c5')),...
    cat(2,getd(evt,'d1'),getd(evt,'d3'),getd(evt,'d5'))};
cfg.twin = [-.5,3.5]; 
% filter LFP at 70 - 90 Hz passband
cfg.filter.type = 'cheby1';
cfg.filter.order = 5;
cfg.filter.f = [70,90];
cfg.filter.verbose = 0;
% threshold LFP after filtering
cfg.t.method = 'zscore';
cfg.t.threshold = 3;
cfg.t.operation = '>';
cfg.t.merge_thr = .05;
cfg.t.minlen = .05;
cfg.t.verbose = 0;
% plot the detected events on top of the LFP trace
cfg.plot = 2; 
eventLFPplot(cfg,csc)







