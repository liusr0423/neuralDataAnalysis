function ftrl_evt = eventLFPfilter(cfg_in,trl)
% function eventLFPfilter filters LFP in certain frequency band and output
% intervals of event detected after filtering
% 
% REQUIRED input: 
%                  trl:   tsd with LFP data
%                  cfg.f: filter range to use
% CFG options with defaults: 
%                  cfg.ftype = 'butter'; {'cheby1','butter','fdesign'}
%                  cfg.threshold = 3;
%                  cfg.operation =  '>';
%                  cfg.merge_thr = 0.05; % merge events closer than this
%                  cfg.minlen = 0.05; % minimum interval length
%
%  OUTPUT:         
%                  ftrl_evt: time intervals of events detected after 
%                            filtering   
% Sirui Liu Feb-22-2016  

if isempty(cfg_in.f)
    error('Please specify a filter range to use (Hz)! Ex.cfg.f = [140 200];')
end

% filter trial LFP in frequency band specified
ftrl = FilterLFP(cfg_in,trl);

% obtain power and z-score it
ftrl_p = LFPpower([],ftrl);
ftrl_z = zscore_tsd(ftrl_p);

% detect events by thresholding
cfg_def.method = 'raw';
cfg_def.threshold = 3;
cfg_def.operation =  '>'; % return intervals where threshold is exceeded
cfg_def.merge_thr = 0.05; % merge events closer than this
cfg_def.minlen = 0.05;    % minimum interval length
mfun = mfilename;
cfg_in = ProcessConfig(cfg_def,cfg_in,mfun);  % take whatever is in cfg_in and put it into cfg!

ftrl_evt = TSDtoIV(cfg_in,ftrl_z);
        
end
