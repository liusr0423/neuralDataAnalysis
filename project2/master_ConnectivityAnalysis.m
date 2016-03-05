%% some initial path settings
restoredefaultpath; % start with a clean slate
clear classes; 

cd('/Users/sirui/Documents/MATLAB/fieldtrip'); % replace paths with yours
ft_defaults;

cd('/Users/sirui/Documents/MATLAB/neuraldata-w16/shared'); % replace paths with yours
p = genpath(pwd); % create list of all folders from here
addpath(p);

set(0,'DefaultAxesFontSize',18)

%% cd to data folder
df = 'R020';
dataDir = ['../../class/neuralDataAnalysis/Data/',df]; % replace data paths with yours
cd(dataDir);
%%
dataFiles = dir([df,'*']);
numSessions = length(dataFiles);

this_session = dataFiles(7).name;
cd(this_session);
LoadExpKeys;
ExpKeys
%% 
cfg = [];
cfg.fc = {'R020-2012-12-17-CSC03a.ncs','R020-2012-12-17-CSC03d.ncs'};
cfg.label = {'3a','3d'};
csc = LoadCSC(cfg);
 
csc = restrict(csc,ExpKeys.TimeOnTrack(1),ExpKeys.TimeOffTrack(2)); % restrict to task
%%
Fs = csc.cfg.hdr{1}.SamplingFrequency; 
wsize = 2048;

nS = length(csc.label);
for iS = 1:nS
    [P{iS},F{iS}] = pwelch(getd(csc,csc.label{iS}),hanning(wsize),wsize/2,2*wsize,Fs);
 
    for iS2 = iS+1:nS
        [C{iS,iS2},Fc{iS}] = mscohere(getd(csc,csc.label{iS}),getd(csc,csc.label{iS2}),hanning(wsize),wsize/2,2*wsize,Fs);
    end
end
%% 
% plot
subplot(121)
cols = 'kgm';
for iS = 1:nS
h(iS) = plot(F{iS},10*log10(P{iS}),cols(iS),'LineWidth',2); hold on;
end
set(gca,'XLim',[0 150],'XTick',0:5:150,'FontSize',12); grid on;
legend(h,csc.label,'Location','Northeast'); legend boxoff;
xlabel('Frequency (Hz)'); ylabel('Power (dB)'); 

subplot(122); clear h;
h(1) = plot(Fc{1},C{1,2},'LineWidth',2); hold on;
set(gca,'XLim',[0 150],'XTick',0:5:150,'FontSize',12); grid on;
legend(h,{'3a-3d'},'Location','Northeast'); legend boxoff;
xlabel('Frequency (Hz)'); ylabel('Coherence');



