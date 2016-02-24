function eventLFPplot(cfg,csc)
% EVENTLFPPLOT plots the event-triggered LFP traces of data csc
% with the options of 
% 1. Decimate:         resample data at 1/R times the original sample 
%                      rate after lowpass filtering before plotting 
% 2. Signal Filtering: filtering the LFP of each trial before plotting
% 3. Event Detection:  thresholding the filtered data and plot the detected 
%                      intervals on top of the LFP traces after 
% INPUTS required:  
%                  csc:            tsd with LFP data
%                  cfg.eventTimes: cell array of timestamps of events of 
%                                  interest
%                  
%        optional: 
%                  cfg.eventLabel: corresponding event labels to put on
%                                  plot
%                  cfg.twin:       time window of interest (seconds 
%                                  relative to the event times), defualt 
%                                  is [-1 3]
%                  cfg.filter.f: filter options with defaults from FilterLFP()
%                  cfg.filter.t: threshold options with defaults modified from TSDtoIV()
%                  cfg.plot: if 1 plot the detected events on top of the LFP trace, defualt = 0
%                  cfg.r: resample data at cfg.r times of the 
%                         original sample rate, if decimate
%                
% 
% Sirui Liu  Feb-22-2016  

if ~isfield(cfg.twin)
    twin = [-1,3];
else 
    twin = cfg.twin;
end

for ii = 1:length(cfg.eventTimes) % for each event
    figure(ii)
    
    % extract trial timestamps of this event
    trials = cfg.eventTimes{ii};
    
    for n = 1:length(trials) % for each trial of this event
        
        % extract the corresponding piece of LFP according to time window
        % specified of this trial
        trl = restrict(csc,trials(n) + twin(1),trials(n) + twin(2)); 
        
        % filtering tril LFP if specified
        if isfield(cfg,'filter')
         
           [trl,ftrl_evt] = eventLFPfilter(cfg.filter,cfg.threshold,trl);
        end
        
        % decimate trial LFP if specified
        if isfield(cfg,'decimate')
            trl.data = decimate(trl.data,cfg.r);
            trl.tvec = trl.tvec(1:length(trl.data));
        end
        
        % replace the original time axis with a new one based on the time
        % window asked for 
        trl.tvec = trl.tvec - (trl.tvec(1) - twin(1));
        
        % range and mean of LFP plotting
        lfp_minmax = 1; 
        lfp_cent = n;
        
        % rescale the LFP so easier to visualize
        trl.data = rescale(trl.data,-lfp_minmax,lfp_minmax); 
        
        % add a y-offset to the LFP to plot one above the other
        trl.data = trl.data + lfp_cent;    
        
        % plot the event-triggered LFP traces
        plot(trl.tvec,trl.data,'k');
        hold on
        
        % add event label to plot, if any
        if ~isempty(cfg.eventLabel)
            title(cfg.eventLabel{ii},'FontSize',20);
        end
        
        % plot detected intervals after filtering on top of the LFP traces
        if isfield(cfg,'plot') && cfg.plot == 1
            
            % replace the original time axis with a new one based on the
            % time window asked for
            ftrl_evt.tstart = ftrl_evt.tstart  - (trl.tvec(1) - twin(1));
            ftrl_evt.tend = ftrl_evt.tend - (trl.tvec(1) - twin(1));
            
            PlotTSDfromIV([],ftrl_evt,trl);  
       
        end

    end
    
    % add a vertical line to the plot indicating time zero
    y1 = get(gca,'ylim');
    plot([0,0],y1)
    set(gca,'FontSize',20);
    axis xy; xlabel('time (s)'); ylabel('trial');

end

