%Example of How to Use the Filter Function

%Parameters

Fs = 100;
tmax = 5;
Nsamps = tmax*Fs;

%Create Initial Signals
t = 1/Fs:1/Fs:tmax;
s1 = 10*cos(2*pi*t);
s2 = 2*cos(20*pi*t + pi/4);
s3 = s1 + s2;

%Plot in Time Domain

%Original
figure
plot(t,s1)
xlabel('Time (s)')
ylabel('Amplitude (V)')
title('Original Signal')
ylim([-15 15])

%Original + High Freq
figure
plot(t,s3)
xlabel('Time (s)')
ylabel('Amplitude (V)')
title('Original Signal Combined With High Frequency Signal')
ylim([-15 15])

%Filter Signals

%Simple Low-Pass Filter
b = 1;
a = [1 -1];

%Apply Filter
s3_f = filter(b,a,s3);

%Scale Output
s3_f = s3_f/15;

%Plot Filtered Signal
figure
plot(t,s3_f)
xlabel('Time (s)')
ylabel('Amplitude (V)')
title('Filtered Signal')
ylim([-15 15])
%%
Nsamps = Fs * twin(2);

%Frequency Domain

f = [-Nsamps/2:Nsamps/2-1]./Nsamps;   %Prepare freq data for plot

%Original + High Freq
s3_fft = abs(fft(x));
s3_fft =fftshift (s3_fft);      

figure
plot(f, s3_fft)
xlabel('Frequency (Fs^{-1})')
ylabel('Amplitude')
title('Frequency Response of Combined Signal Before Filtering')
ylim([0 100])

%Filter Signals
s3_f_fft = abs(fft(xf));
s3_f_fft = fftshift(s3_f_fft)    ; 

figure
plot(f, s3_f_fft)
xlabel('Frequency (Fs^{-1})')
ylabel('Amplitude')
title('Frequency Response of Combined Signal After Filtering')
ylim([0 100])