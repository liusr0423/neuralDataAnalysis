%% set MATLAB's current directory to the data folder (R042-2013-08-18);
fd = '/Users/liusirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R042-2013-08-18';
cd(fd);
%% load the data (note, may need to unzip position data first)
cfg = [];
S = LoadSpikes(cfg);

%% load a LFP
cfg = [];
cfg.fc = {'R042-2013-08-18-CSC03a.ncs'};
csc = LoadCSC(cfg); % continuously sampled channel

%% restrict the data to interval of interest
this_iv = iv([5900 6000]);
S_r = restrict(S,this_iv);
csc_r = restrict(csc,this_iv);

%% plot spikes
SET_spkY = 0.4; % parameter to set spike height
figure; hold on;
for iC = 1:length(S_r.t)
    if ~isempty(S_r.t{iC})
        plot([S_r.t{iC} S_r.t{iC}],[iC-SET_spkY iC+SET_spkY],'k');
    end
end % of cells
ylabel('neuron #');