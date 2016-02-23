function ftrl_evt = eventLFPfilter(cfg_f,cfg_t,trl)
% EVENTLFPFILTER filters LFP in certain frequency band and output
% time intervals of intervals detected after filtering
% 
% REQUIRED input: 
%                  trl:   tsd with LFP data
%                  cfg_f.f: filter range to use
% CFG_f filter options with defaults from function FilterLFP():
%                  cfg_f.order = 4; % filter order
%                  cfg_f.display_filter = 0; % show output of fvtool on filter
%                  cfg_f.bandtype = 'bandpass'; % 'highpass', 'lowpass'
%                  cfg_f.R = 0.5; % passband ripple (in dB) for Chebyshev filters only
%                  cfg_f.ftype = 'butter'; {'cheby1','butter','fdesign'}
%                  cfg_f.verbose = 1; If 1 display helpful text in command window, if 0 don't
% CFG_t threshold options with defaults modified from TSDtoIV(): 
%                  cfg_t.method = 'raw';
%                  cfg_t.threshold = 3;
%                  cfg_t.operation =  '>';
%                  cfg_t.merge_thr = 0.05; % merge events closer than this
%                  cfg_t.minlen = 0.05; % minimum interval length
%                  cfg.minlen = .05;
%                  cfg.verbose = 1; 1 display command window, 0 don't
%
%  OUTPUT:         
%                  ftrl_evt: time intervals of events detected after 
%                            filtering   
% 
% Feb-22-2016 by Sirui Liu 

if isempty(cfg_f.f)
    error('Please specify a filter range to use (Hz)! Ex.cfg.f = [140 200];')
end

% filter trial LFP in frequency band specified
ftrl = FilterLFP(cfg_f,trl);

% obtain power and z-score it
ftrl_p = LFPpower([],ftrl);
ftrl_z = zscore_tsd(ftrl_p);

% detect events by thresholding
cfg_def.method = 'raw';
cfg_def.threshold = 3;
cfg_def.operation =  '>'; % return intervals where threshold is exceeded
mfun = mfilename;
cfg_t = ProcessConfig(cfg_def,cfg_t,mfun);  % replace default cfg_def with cfg_t
ftrl_evt = TSDtoIV(cfg_t,ftrl_z);
        
end
