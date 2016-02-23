function eventLFPplot(cfg,csc)
% function eventLFPplot plots the event-triggered LFP traces of data csc
% INPUTS required:  
%                  csc: tsd with LFP data
%                  cfg.eventTimes: cell array of timestamps of events of 
%                                  interest
%                  
%        optional: cfg.twin: time window of interest (seconds relative to 
%                            the event times), defualt is [-1 3
%                  cfg.f: filter range to use, if filtering LFP to detect events 
% 
% Sirui Liu  Feb-22-2016  

if isempty(cfg.twin)
    twin = [-1,3];
else 
    twin = cfg.twin;
end


for ii = 1:length(cfg.eventTimes) % for each event
    figure
    
    % extract trial timestamps of this event
    trials = cfg.eventTimes{ii};
   
    for n = 1:length(trials) % for each trial of this event
        
        % extract the corresponding piece of LFP according to time window
        % specified of this trial
        trl = restrict(csc,trials(n) + twin(1),trials(n) + twin(2)); 
        
        % filter trial LFP in frequency band as specified in cfg.f, if any
        if cfg.f
            ftrl_evt = eventLFPfilter(struct('f',cfg.f),trl);
            
            % replace the original time axis with a new one based on the
            %  time window asked for
            ftrl_evt.tstart = ftrl_evt.tstart  - (trl.tvec(1) - twin(1));
            ftrl_evt.tend = ftrl_evt.tend - (trl.tvec(1) - twin(1));
        end
        
        % replace the original time axis with a new one based on the time 
        % window asked for 
        trl.tvec = trl.tvec - (trl.tvec(1) - twin(1));
        
        % range and mean of LFP plotting
        lfp_minmax = 1; lfp_cent = n;
        
        % rescale the LFP so easier to visualize
        trl.data = rescale(trl.data,-lfp_minmax,lfp_minmax); 
        
        % add a y-offset to the LFP to plot one above the other
        trl.data = trl.data + lfp_cent;    
        
        % plot the event-triggered LFP traces
        plot(trl.tvec,trl.data,'k');
        hold on
        if cfg.f
            PlotTSDfromIV([],ftrl_evt,trl);    % plot detected events after filtering on top of the lfp plot
        end
    end
    
    % add a vertical line to the plot indicating time zero
    y1 = get(gca,'ylim');
    plot([0,0],y1)
    set(gca,'FontSize',20);
    axis xy; xlabel('time (s)'); ylabel('trial');

end

