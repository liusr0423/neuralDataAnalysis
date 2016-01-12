%% load data
cd('/Users/liusirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R016/R016-2012-10-08_promoted'); % same session as Module 1
 
cfg = [];
cfg.fc = {'R016-2012-10-08-CSC02d.ncs'}; % cell array with filenames to load
csc = LoadCSC(cfg);

%% remember to use Cell Mode in the editor to run this code! 
cfg = [];
evt = LoadEvents(cfg);

%% load video data (make sure the VT1.zip file is unzipped first and now present in MATLAB's working folder!)
[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, [] );

%% plot video data -- use a new cell so that you can rerun this without also reloading the data
fh = figure; set(fh,'Color',[0 0 0]);
keep_idx = X; keep_idx(keep_idx==0)=[];
keep_idy = Y;  keep_idy(keep_idy==0)=[];
plot(X(keep_idx),Y(keep_idx),'.','Color',[0.7 0.7 0.7],'MarkerSize',1); axis off;