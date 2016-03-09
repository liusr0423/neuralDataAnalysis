cfg              = [];
cfg.output       = 'powandcsd';  % output power and cross-spectral density
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100;      % frequencies to use
cfg.t_ftimwin    = ones(length(cfg.foi));  % frequency-dependent, 20 cycles per time window
cfg.keeptrials   = 'yes';
cfg.channel      = {'HC' 'Str'};
cfg.channelcmb   = {'HC' 'Str'}; % channel pairs to compute csd for
cfg.toi          = -2:0.05:4.5;
TFR              = ft_freqanalysis(cfg, data_trl);

cfg_psi = [];
cfg_psi.method = 'psi';
cfg_psi.bandwidth = 4; % number of frequencies to compute slope over
cfg_psi.channel = {'HC','Str'};
cfg_psi.channelcmb = {'Str','HC'; 'HC','Str' };
C = ft_connectivityanalysis(cfg_psi,TFR);

% plot
figure
lbl = [C.labelcmb{1,:}]; % get the label of this pair
imagesc(C.time,C.freq,sq(C.psispctrm(1,:,:))); axis xy; colorbar; colormap jet
h = colorbar;
ylabel(h, 'psi(?/Hz)')
xlabel('time (s)'); ylabel('Frequency (Hz)'); title(lbl);