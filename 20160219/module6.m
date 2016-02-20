%% using filter()
load count.dat;
x = count(:,1);
a = 1; % a_0 is the (hidden) coefficient on the left side, in front of y(n)
b = [1/4 1/4 1/4 1/4]; % four b's of 1/4 each so we get the mean
y = filter(b,a,x); % x is the original signal, y the filtered version

t = 1:length(x);
plot(t,x,'-.',t,y,'-'), grid on
legend('Original','Filtered',2)
% from the plot, first sample equals x(1)/4, input outside assumes to be 0 (edge effect)
% the filtered signal is phase-shifted toward the right this of course 
% arises because our y(n) is based only on past samples, not on the future.
%% exercise using Butterworth filter
% set up time axis
Fs = 500; %Hz
dt = 1./Fs;
twin = [0 10]; 
tvec = twin(1):dt:twin(2)-dt; 

%generate white noise
x = rand(1,length(tvec));

% get original PSD in dB using a 512-samping Hanning window evaluated over
% 2^14 points
wSize = 512;
nP = 2^14;
[Porig,Forig] = pwelch(x,hanning(wSize),wSize/2,nP,Fs);

% design filter
W1 = 50; %Hz
W2 = 100;
W1 = W1 / (Fs/2); % normalize by Fs/2
W2 = W2 / (Fs/2);
order = 4;
[b,a] = butter(order,[W1,W2]); %bandpass Butterworth filter in order of 4
y = filter(b,a,x);
% get filtered PSD
[Pfilt,Ffilt] = pwelch(y,hanning(wSize),wSize/2,nP,Fs);

figure
set(gcf,'color','w')
subplot(1,2,1)
plot(Forig,10*log10(Porig),'b'); 
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
set(gca,'YLim',[-140,-20])
set(gca,'XLim',[0,250])
set(gca,'XTick',0:50:250)
title('original')
grid on
subplot(1,2,2)
plot(Ffilt,10*log10(Pfilt),'b'); 
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
set(gca,'YLim',[-140,-20])
set(gca,'XLim',[0,250])
set(gca,'XTick',0:50:250)
title('filtered')
grid on

%% use buttord() to find an appropriate filter order
% The buttord() function takes the filter specifications and returns a 
% suggested filter order (N) and new frequency cutoffs Wn to feed to butter(). 
Wp = [ 50 100] * 2 / Fs; % passband - between 50 and 100 Hz
Ws = [ 45 105] * 2 / Fs; % stopband
Rp = 3;  % passband ripple of no more than Rp
Rs = 20; % minimum level of attenuation in the stopband
[N,Wn] = buttord( Wp, Ws, Rp, Rs); % determine filter parameters
[b2,a2] = butter(N,Wn); % builds filter
fvtool(b,a,b2,a2); % Filter Visualization Tool 
legend('old filter','new filter')

%% Chebyshev Type I Filter

Wp = [ 50 100] * 2 / Fs; 
Ws = [ 48 102] * 2 / Fs;
[N,Wn] = cheb1ord( Wp, Ws, 3, 20); 
[b_c1,a_c1] = cheby1(N,0.5,Wn);
H = fvtool(b2,a2,b_c1,a_c1);
legend('Butterworth filter','Chebyshev filter')
zoom(H, [.15,.45,-35 0]) %constrain the zoom to the y-axis.
% sharper but passband is not flat


%% phase responses 
% phase shift 
Fs = 500; dt = 1./Fs;
t = [0 10];
tvec = t(1):dt:t(2)-dt;
f1 = 80;
f2 = 40;
s1 = sin(2*pi*f1*tvec+pi/6);
s2 = sin(2*pi*f2*tvec);
s = s1 + s2;
 
sf = filter(b_c1,a_c1,s); % highpass filter(50-100)
figure 
set(gcf,'color','w')
plot(tvec,s,'k',tvec,sf,'r--'); hold on;
legend({'original','filtered'});
xlim([0 0.2]);

%% fillfilt 
% filter the singal forwards and backwards so net phase response is
% zero

sf = filtfilt(b_c1,a_c1,s);
figure
plot(tvec,s,'k',tvec,sf,'r--'); hold on;
legend({'original','filtered'});
xlim([0 0.2]);
 
%% compare freq responses
% original vs. chebyshev filter vs. filtfiilt 
Fs = 500; dt = 1./Fs;
t = [0 10];
tvec = t(1):dt:t(2)-dt;
 
x = rand(size(tvec)); % white noise input
[P,F] = pwelch(x,hanning(512),256,2^14,Fs);
 
y1 = filter(b_c1,a_c1,x);
[P1,F1] = pwelch(y1,hanning(512),256,2^14,Fs);
 
y2 = filtfilt(b_c1,a_c1,x);
[P2,F2] = pwelch(y2,hanning(512),256,2^14,Fs);

figure
set(gcf,'color','w')
plot(F,10*log10(P),F,10*log10(P1),F,10*log10(P2));
legend({'original','filter','filtfilt'});
% filtfilt is steeper since filtering twice (~increase the order of the
% filter)

%% exercise
% find out if decimate() produces phase shifts
Fs = 500; % Hz
t = 0:1/Fs:1-1/Fs;  % Time vector
x = sin(2*pi*30*t) + sin(2*pi*60*t);
y_downsample = downsample(x,4);
y = decimate(x,4); % decimate function applies low-pass chebyshev with 
                   % cut-off freq .8*(Fs/2)/R before resampling
figure
set(gcf,'color','w')
plot(y_downsample,'b:')
axis([0 125 -2 2] )% Original signal desampled
hold on;
plot(y,'r-')                        % Decimated signal
legend('Original Signal','Decimated Signal')

%% neuroscience applications
% notch filter to remove 60hz line noise
%% old method
Fs = 500; %Hz
[b,a] = butter(10, [59 61] * 2 / Fs, 'stop');
H = fvtool(b,a);
zoom(H,'x',[.2,.3])

%% new method
Fs = 500; %Hz
[z,p,k] = butter(10, [59 61] * 2 / Fs, 'stop'); % note, we ask for 3 outputs instead of 2
[sos,g] = zp2sos(z,p,k); % convert to SOS format
h = dfilt.df2sos(sos,g); % create filter object: 'second-order section' format 
fvtool(h);

%% test on white noise using filtfilt()
%generate white noise
Fs = 500; %Hz
dt = 1./Fs;
twin = [0 10]; 
tvec = twin(1):dt:twin(2)-dt; 
x = rand(size(tvec));
% filter with notch filter to remove 60hz line noise
xf = filtfilt(sos,g,x);

figure 
set(gcf,'color','w')
plot(tvec,x,'k',tvec,xf,'r-'); hold on;
legend({'original','filtered(band-stop 60hz)'});
xlim([0 0.4]);
%% 
% get original PSD in dB using a 512-samping Hanning window evaluated over
% 2^14 points
wSize = 512;
nP = 2^14;
[Porig,Forig] = pwelch(x,hanning(wSize),wSize/2,nP,Fs);

% get filtered PSD
[Pfilt,Ffilt] = pwelch(xf,hanning(wSize),wSize/2,nP,Fs);

figure
set(gcf,'color','w')
subplot(1,2,1)
plot(Forig,10*log10(Porig),'b'); 
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
set(gca,'YLim',[-140,-20])
set(gca,'XLim',[0,250])
set(gca,'XTick',0:50:250)
title('original')
grid on
subplot(1,2,2)
plot(Ffilt,10*log10(Pfilt),'b'); 
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
set(gca,'YLim',[-140,-20])
set(gca,'XLim',[0,250])
set(gca,'XTick',0:50:250)
title('filtered')
grid on

%% test on another signal
Fs = 500; dt = 1./Fs;
t = [0 10];
tvec = t(1):dt:t(2)-dt;
f1 = 60;
f2 = 30;
s1 = sin(2*pi*f1*tvec+pi/6);
s2 = sin(2*pi*f2*tvec);
s = s1 + s2;

sf = filtfilt(sos,g,s);


figure 
set(gcf,'color','w')
plot(tvec,s,'k',tvec,sf,'r-'); hold on;
legend({'original','filtered(band-stop 60hz)'});
xlim([0 0.4]);

%% Detecting movement artifacts
% cd to your location here
rootDir = '/Users/sirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R016/R016-2012-10-08_promoted';
cd(rootDir);
cfg = [];
cfg.fc = {'R016-2012-10-08-CSC02b.ncs'};
csc = LoadCSC(cfg);
 
cscR = restrict(csc,1270,1272);
plot(cscR.tvec,cscR.data)


%% filtering chewing events
Fs = cscR.cfg.hdr{1}.SamplingFrequency;
Wp = [ 180 220] * 2 / Fs; % pass band
Ws = [ 178 222] * 2 / Fs; % stop band
[N,Wn] = cheb1ord( Wp, Ws, 3, 20); % determine filter parameters
[b_c1,a_c1] = cheby1(N,0.5,Wn); % builds filter
 
fvtool(b_c1,a_c1); % remember to check your filter!
%%  pick up chewing events
y = filtfilt(b_c1,a_c1,cscR.data);
plot(cscR.tvec,cscR.data,'b',cscR.tvec,y,'r');
chew_power = y.^2; %signal power
chew_power_filtered = medfilt1(chew_power,101); % filter window is specified in samples, so this is ~50ms(1/Fs*101)
[h1 h2] = plotyy(cscR.tvec,cscR.data,cscR.tvec,chew_power_filtered);


%% detect SWR events

rootDir = '/Users/sirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R042-2013-08-18';
cd(rootDir);
cfg = []; cfg.fc = {'R042-2013-08-18-CSC03a.ncs'};
lfp = LoadCSC(cfg);
%% filter in SWR band
cfg = [];
cfg.f = [140 220]; % filter range to use (Hz)
cfg.display_filter = 0; % show output of fvtool on filter
SWRf = FilterLFP(cfg,lfp); %input: tsd with LFP data; output: tsd with filtered data

%% obtain power and z-score it
SWRp = LFPpower([],SWRf); % input: tsd with filtered LFP data; output: tsd with power envelope obtained by Hilbert transform
SWRp_z = zscore_tsd(SWRp); %   z-scores data

%% detect events
cfg = [];
cfg.method = 'raw';
cfg.threshold  = 3;
cfg.operation = '>'; % return intervals where threshold is 
cfg.merge_thr = .05; % merge event closer than this
cfg.minlen = .05; % minimum interval length
%  create interval data from tsd by tresholding
SWR_evt = TSDtoIV(cfg,SWRp_z);

%% to each event, add a field with the max z-scored power (for later selection)
cfg = [];
cfg.method = 'max'; % 'min', 'mean'
cfg.label = 'maxSWRp'; % what to call this in iv, i.e. usr.label
% AddTSDtoIV(cfg_in,iv_in,tsd_in): add usr field to iv based on tsd
SWR_evt = AddTSDtoIV(cfg,SWR_evt,SWRp_z); 

%% select only those events of >5 z-scored power
cfg = [];
cfg.operation = '>';
cfg.threshold = 5;
SWR_evt = SelectIV(cfg,SWR_evt,'maxSWRp'); % select iv from iv
%% plot events in highlighted on top of full lfp
PlotTSDfromIV([],SWR_evt,lfp); % display intervals defined relative to tsd
%% ..or the events alone (fixed 200ms window) centered at event time)
close all;
cfg = [];
cfg.display = 'iv';%only plot iv's
cfg.mode = 'center'; % plot each iv using specific time window ('center')
cfg.fgcol = 'k';% tsd color outside iv's
PlotTSDfromIV(cfg,SWR_evt,lfp);
%% ..hold on (highlight edges of event on top of previous plot)
cfg = [];
cfg.display = 'iv';%only plot iv's
cfg.fgcol = 'r';%tsd color within iv's
PlotTSDfromIV(cfg,SWR_evt,lfp);