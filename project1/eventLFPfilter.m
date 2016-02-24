function [ftrl,ftrl_evt] = eventLFPfilter(cfg_f,cfg_t,trl)
% EVENTLFPFILTER filters LFP and threshold the filtered data (optional)  
% 
% REQUIRED input: 
%                  trl:     tsd with LFP data
% CFG_f filter options with defaults from function FilterLFP():
%                  cfg_f.order = 4; % filter order
%                  cfg_f.display_filter = 0; % show output of fvtool on filter
%                  cfg_f.bandtype = 'bandpass'; % 'highpass', 'lowpass'
%                  cfg_f.f = [6 10]; filter range to use (in Hz)
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
%                  ftrl : filtered LFP
%                  ftrl_evt: time intervals of events detected after 
%                            filtering   
% 
% Feb-22-2016 by Sirui Liu 

% filter trial LFP
cfg_f.verbose = 0;
ftrl = FilterLFP(cfg_f,trl);

if cfg_t % if thresholding for event detection
    
    % obtain power and z-score it
    ftrl_p = LFPpower([],ftrl);
    ftrl_z = zscore_tsd(ftrl_p);
    ftrl_evt = TSDtoIV(cfg_t,ftrl_z);

end

end
