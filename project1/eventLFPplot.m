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
%                  cfg.filter.f:   filter configurations with defaults from
%                                  FilterLFP()
%                  cfg.filter.t:   threshold configurations with defaults 
%                                  from TSDtoIV()
%                  cfg.plot:       default = 0: plot the original LFP trace
%                                  if 1 plot the filtered LFP trace
%                                  if 2 plot the detected events on top of
%                                  the original LFP trace
%                                  if 3 plot the detected events on top of
%                                  the filtered LFP trace
%                  cfg.r:          resample data at cfg.r times of the
%                                  original sample rate, if decimate
%
%
% Sirui Liu  Feb-22-2016

if ~isfield(cfg,'twin')
    twin = [-1,3];
else
    twin = cfg.twin;
end

if ~isfield(cfg,'plot')
    cfg.plot = 0;
end

for ii = 1:length(cfg.eventTimes) % for each event
    
    % extract trial timestamps of this event
    trials = cfg.eventTimes{ii};
    
    for n = 1:length(trials) % for each trial of this event
        
        % extract the corresponding piece of LFP according to time window
        % specified of this trial
        trl = restrict(csc,trials(n) + twin(1),trials(n) + twin(2));
        
        % filtering trial LFP if specified
        if isfield(cfg,'filter')
           
            ftrl = FilterLFP(cfg.filter,trl); % filter

            if isfield(cfg,'t') % thresholding for event detection
                
                ftrl_p = LFPpower([],ftrl);       % obtain power envelope
                ftrl_evt = TSDtoIV(cfg.t,ftrl_p); % detect events

                % replace the original time axis of detected events with a
                % new one based on the time window asked for
                ftrl_evt.tstart = ftrl_evt.tstart - (ftrl.tvec(1) - ...
                    twin(1));
                ftrl_evt.tend = ftrl_evt.tend - (ftrl.tvec(1) - twin(1));
            end
        end
        
        if  isfield(cfg,'plot') && mod(cfg.plot,2)
            trl = ftrl;
        end
        
        % decimate trial LFP if specified
        if isfield(cfg,'r')
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
        
        % plot only the event-triggered LFP traces
        figure(ii)
        plot(trl.tvec,trl.data,'k');
        hold on
        % plot detected intervals on top of the LFP traces
        if cfg.plot == 2 || cfg.plot == 3
            PlotTSDfromIV([],ftrl_evt,trl);
        end
        
        % add event label to plot, if any
        if ~isempty(cfg.eventLabel)
            title(cfg.eventLabel{ii},'FontSize',20);
        end
        
    end
    
    % add a vertical line to the plot indicating time zero
    y1 = get(gca,'ylim');
    plot([0,0],y1)
    set(gca,'FontSize',20);
    xlabel('time (s)'); ylabel('trial');
    set(gca,'XTick',[twin(1):.5:twin(2)]);
    axis xy;
    
end

