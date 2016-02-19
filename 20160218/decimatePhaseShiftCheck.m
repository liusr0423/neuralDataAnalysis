%% downsampling
s_100_1 =  makesin(100,1,1);
s_20_1 = makesin(20,1,1);
x = s_100_1(1:5:100);% just downsampling by picking every 5th sample...
max(abs(x - s_20_1)) % works fine!
subplot(2,1,1)
plot(x,'b-')
hold on
plot(s_20_1,'r:')
title('downsample')
legend('downsampled from 100hz','sine wave(20hz)')

%% decimate
x = decimate(s_100_1,5);
max(abs(x - s_20_1))
subplot(2,1,2)
plot(x,'b-')
hold on
plot(s_20_1,'r:')
title('decimate')
legend('decimated from 100hz','sine wave(20hz)')
