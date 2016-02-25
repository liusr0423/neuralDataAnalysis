#Project 1: Spectrogram (module 7)

Codes used for plotting averaged event-triggered spectrogram and corresponding trial LFP traces for each event type defined on data session R016-2012-10-30-CSC04a from vandermeer lab server.

Run master_Spectrogram.m (replace code and data paths with yours) with following steps:
- Some initial path settings
- Load data files from session R016-2012-10-03-CSC04a
- Load events info for this session
- For each event type, plot the averaged event-triggered spectrogram (frequency-dependent windowing) over all trials
- For each event type, plot the event-triggered trial LFP traces using eventLFPplot(cfg,csc) with the options of 
   1.  Decimate
   2.  Signal Filtering
   3.  Event Detection

Used MATLAB R2015a running on Mac OS X
