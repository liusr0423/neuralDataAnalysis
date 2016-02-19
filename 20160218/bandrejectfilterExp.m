clear all
alphap=2;  %passband attenuation in dB
alphas=20; %stopband attenuation in dB
wp=[.2*pi,.4*pi]; % passband freq. in radians
ws=[.1*pi,.5*pi]; % stopband freq. in radians
%to find cutoff freq. and order of the filter
[n,wn]=buttord(wp/pi,ws/pi,alphap,alphas); %syatem function of the filter
[b,a]=butter(n,wn,'STOP');
w=0:.01:pi;

[h,ph]=freqz(b,a,w);
m=20*log(abs(h));
an=angle(h);
subplot(2,1,1);plot(ph/pi,m);grid;
ylabel('Gain in dB');
xlabel('NORMALISED FREQUENCY')
subplot(2,1,2);plot(ph/pi,an);grid
ylabel('phase in radians');
xlabel('noramlised frequency')
