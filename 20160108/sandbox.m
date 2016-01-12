%% load data
cd('/Users/sirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R016/R016-2012-10-08_promoted'); % replace this with where you saved the data
 
cfg = [];
cfg.fc = {'R016-2012-10-08-CSC02d.ncs'}; % cell array with filenames to load
csc = LoadCSC(cfg);

%%
plot(csc.tvec,csc.data);
xlim([1338.6 1339.2]);