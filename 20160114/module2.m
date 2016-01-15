%% Timestamped data (TSD) data-type: load data
cd('/Users/liusirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R016/R016-2012-10-08_promoted'); % same session as Module 1
 
cfg = [];
cfg.fc = {'R016-2012-10-08-CSC02d.ncs'}; % cell array with filenames to load: 
                                         % a single local field potential, recorded 
                                         % from a specific electrode in the brain
csc = LoadCSC(cfg);

%% plot this data
plot(csc.tvec,csc.data);
%% how does LoadCSC() represent multiple signals
cfg.fc = FindFiles('*CSC01*.ncs');
csc = LoadCSC(cfg);
size(csc.data) %  nSignals x nSamples 
%%  plot out specific channel using getd on tsd data
plot(getd(csc,'R016-2012-10-08-CSC01a.ncs'),0,'.k') 

%% Timestamp (TS) data-type: remember to use Cell Mode in the editor to run this code! 
cfg = [];
evt = LoadEvents(cfg);
%% retrieve times associated with 1 pellet cue and plot each time against zero
plot(getd(evt,'1 pellet cue'),0,'.k') 
%% how to use LoadSpikes to represent multiple ts neurons
S = LoadSpikes([]); % This instructs LoadSpikes() to load all spike files it can find
%% Interval (IV) data-type: 
a = iv([1 2]); % define a single interval from 1 to 2
b = iv([1 2],[3 3]); % define two intervals, 1 to 3 and 2 to 3
