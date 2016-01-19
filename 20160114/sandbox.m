%% set MATLAB's current directory to the data folder (R042-2013-08-18);
cd('/Users/liusirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R042-2013-08-18'); 

%% load video data 
% (make sure the VT1.zip file is unzipped first and now present in MATLAB's working folder!)
% use low-level loading function Nlx2MatVT for video data
[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, [] );
whos
%% plot it
plot(X,Y);

%% plot it again without the missing data
fh = figure; set(fh,'Color',[0 0 0]); % set background to black
keep_idx = find(X ~=0 | Y ~=0);
plot(X(keep_idx),Y(keep_idx),'.','Color',[0.7 0.7 0.7],'MarkerSize',1); axis off; % set poitns to size 1 in a gray color

%% save figures
set(gcf, 'InvertHardCopy', 'off'); % preserve the black background
print(gcf,'-r75','-dpng','module2_xvsy2.png'); % saves a 75dpi PNG image

%% Timestamps 
plot(Timestamps(keep_idx),X(keep_idx),'.r','MarkerSize',3)
box off;
set(gca,'FontSize',24);

%% LFP data file (*.Ncs) loading
clear all;
fname = 'R042-2013-08-18-CSC05a.ncs';
[Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);
% Neuralynx Ncs data is stored in blocks of 512 samples
% with only the first sample of each block timestamped. 
% Hence the [512 x 17193] size of Samples

%% Event file (*.Nev) loading
fn = FindFile('*Events.nev');
[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);

%% use the wrapped data loaders
% LoadPos()
posdata = LoadPos([]); % note empty config

% plot without knowing which variable is which 
plot(getd(posdata,'x'),getd(posdata,'y'),'.');

%% load csc
cfg = []; % starting with an empty config is good practice -- that way you avoid carryover of previous values!
cfg.fc = {'R042-2013-08-18-CSC05a.ncs'};
csc = LoadCSC(cfg);

%% load events
evt = LoadEvents([]);
%% using cfg to load events
cfg = [];
cfg.eventList = {'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0004).','TTL Output on AcqSystem1_0 board 0 port 0 value (0x0040).'};
cfg.eventLabel = {'FoodDelivery','WaterDelivery'};
evt = LoadEvents(cfg);
%% load spikes
S = LoadSpikes([]);

%% putting it all together: example 1: use of restrict()
LoadMetadata;
 
pos = LoadPos([]);
 
left_pos = restrict(pos,metadata.taskvars.trial_iv_L); % left trials only
plot(getd(left_pos,'x'),getd(left_pos,'y'),'.'); % looks like right trials! camera reverses image, can fix with set(gca,'YDir','reverse')
 
%% example 2: interplay between tsd and iv data
LoadExpKeys;
 
please = []; please.fc = ExpKeys.goodSWR(1); % local field potential with good "sharp wave-ripple" events
lfp = LoadCSC(please); % aacarey is Canadian and asks nicely; cfg name is just arbitrary
 
% detect possible artifacts
cfg = [];
cfg.method = 'zscore'; % first normalize the data
cfg.threshold = -8;
cfg.minlen = 0; % no minimum length on events to detect
cfg.dcn = '<'; % detect intervals with z-score lower than threshold
artifact_iv = TSDtoIV(cfg,lfp); % creates iv with start and end times of possible artifacts
 
% plot detected intervals
cfg = []; cfg.display = 'tsd'; % also try 'iv' mode!
PlotTSDfromIV(lfp,artifact_iv,lfp)