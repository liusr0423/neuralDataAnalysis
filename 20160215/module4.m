%% Detailed examination of Neuralynx time series data
% cd to your location here
rootDir = '/Users/sirui/Documents/MATLAB/class/neuralDataAnalysis/Data/R016/R016-2012-10-08_promoted';
cd(rootDir);
fname = 'R016-2012-10-08-CSC03b.ncs';
% Neuralynx loader
[Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);


%% convert Timestamps to secs
Timestamps_Secs = Timestamps/(10^6);
Samples_num = size(Samples,1) * size(Samples,2);
% the total time that would be sampled continuously if all samples were acquired without any gaps
Time_total = Samples_num / SampleFrequencies(1); % sec
% acutal time elapsed 
Time_total_real = Timestamps_Secs(end) - Timestamps_Secs(1);
%  there are several gaps in the data: plotting the difference between each sample and its predecessor 
plot(diff(Timestamps))

%% create Value session timestamps
LoadExpKeys
ExpKeys
[cON, indexON] = min(abs(Timestamps_Secs - ExpKeys.TimeOnTrack(1)));
[cOFF, indexOFF] = min(abs(Timestamps_Secs - ExpKeys.TimeOffTrack(1)));
TimestampsValue = Timestamps_Secs(indexON:indexOFF);
SamplesValue = Samples(indexON:indexOFF);
NumberOfValidSamplesValue = NumberOfValidSamples(indexON:indexOFF);
plot(diff(TimestampsValue)) % plot and find the most common value between timestamps

%% expected difference between 512-sample timestamps if Fs is 2kHz
512.*(1/2000) == mode(diff(TimestampsValue));
SamplingFreq_true = 512/mode(diff(TimestampsValue))

%% 
plot(diff(TimestampsValue))
hold on;
plot(NumberOfValidSamplesValue == 512,'r')
%the odd timestamp diffs occur for those sample blocks that have a number of valid samples that is not 512.
oddBlocksValue_idx = find(NumberOfValidSamplesValue ~= 512);
oddBlocks_idx = find(NumberOfValidSamples ~= 512);
