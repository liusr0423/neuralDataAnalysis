%% set MATLAB's current directory to the data folder (R042-2013-08-18);
fd = '/Users/liusirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R042-2013-08-18';
cd(fd);
%% load the data (note, may need to unzip position data first)
cfg = [];
S = LoadSpikes(cfg); % S.t contains the spike times from a putative neuron

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
 
figure('Color','white'); hold on;
for iC = 1:length(S_r.t)% plot the spike times for each neuron
    if ~isempty(S_r.t{iC})
        plot([S_r.t{iC} S_r.t{iC}],[iC-SET_spkY iC+SET_spkY],'k');
    end
end % of cells
ylabel('neuron #');
 
%% add a LFP
SET_cscY = [-5 0];
plot(csc_r.tvec,rescale(csc_r.data,SET_cscY(1),SET_cscY(2)),'r');
set(gca,'YTickLabel',{});
 
%% add multi-unit activity in separate axes
ax1 = gca; ax1_pos = get(ax1,'Position');
ax2 = axes('Position',ax1_pos); % new axes with same position as first
 
cfg = []; cfg.tvec = csc.tvec; cfg.sigma = 0.1;
mua = getMUA(cfg,S); % obtain multi-unit activity
 
xr = get(ax1,'XLim');
mua_r = restrict(mua,xr(1),xr(2)); % only keep the data we need
 
axes(ax2); % set current axes to the second set, so what follows happens there
mua_hdl = plot(mua_r.tvec,mua_r.data,'Color',[0.7 0.7 0.7]);
 
set(gca,'YAxisLocation','right','Box','off','XTick',[],'Color','none','YColor',[0.5 0.5 0.5])
%set the properties of the new multi-unit axes to have the y-axis on the right
ylabel('multi-unit activity (spk/s)');
linkaxes([ax1 ax2],'x'); % link the x-axis of both axes so that both update when zooming in
%% zoom in
xlim([5965,5969]);
%You should now be looking at the synchronous activation of a substantial 
%number of neurons, reflected in the MUA peak, and associated with a high-frequency
%oscillation (~150-250Hz) in the LFP. These are the neurophysiological signatures of
%a ?sharp wave-ripple complex? (SWR for short), events which are thought to 
%contribute to the consolidation and retrieval of episodic memories.
%% use the left and right arrow keys to scroll through the rasterplot.
set(gcf,'KeyPressFcn',@figscroll) 

%% movie version
h = gcf; set(h,'Position',[100 100 640 480]);
% keep the size of the resulting movie file manageable 
% (the above sets a 640×480 pixel figure size)
t = [5900 6000]; % start and end times (experiment time)
FPS = 30; % frame rate (per s)
twin = [-1 1]; % width of time window (in s)
tvec = t(1):1/FPS:t(2);
for iT = 1:length(tvec)
   set(gca,'XLim',twin+tvec(iT));
   f(iT) = getframe(gcf); % store current frame
   drawnow; pause(1/FPS);
end


%% using existing visualization tools
%%  multiraster: load data
S = LoadSpikes([]);
please = []; please.fc = {'R042-2013-08-18-CSC03a.ncs'};
csc = LoadCSC(please);
 
LoadMetadata; % load some experiment metadata such as trial start and end times
 
% plot data
cfg_plot = [];
cfg_plot.lfp = csc;
cfg_plot.spkColor = 'jet';
cfg_plot.evt = metadata.taskvars.trial_iv_L; % "left" trials on the T-maze
 
h = MultiRaster(cfg_plot,S);

%% PlotTSDFromIV
cfg = [];
cfg.display = 'iv';
cfg.width = 0.1;
PlotTSDfromIV(cfg,metadata.SWRtimes,csc); % need to LoadMetadata to make this work!

%% fieldtrip
% convert to ft format
cfg = []; cfg.mode = 'resample';
csc_ft = TSDtoFT(cfg,csc);
t0 = csc_ft.time{1}(1); csc_ft.time{1} = csc_ft.time{1}-t0;
 %%
% create ft trials from iv
trl_cfg = [];
trl_cfg.t = IVcenters(metadata.SWRtimes)-t0;
trl_cfg.mode = 'neuralynx';
trl_cfg.hdr = csc_ft.hdr;
trl_cfg.twin = [-1 1];
 
trl = ft_maketrl(trl_cfg);
 
% use trials to create trialified data structure
temp_cfg = []; temp_cfg.trl = trl;
 
ft_in = ft_redefinetrial(temp_cfg,csc_ft);