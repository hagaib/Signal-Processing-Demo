clear
disp('Signal 1 Anti-Aliasing and Resampling');

%% Load Signal Data
outputdir = 'output_data';
filename = 'signal1.mat';
filepath = fullfile(fileparts(mfilename('fullpath')) ,outputdir);
filenamepath = fullfile(outputdir , filename);
load(filenamepath)
obj.plot(1)

%% Signal 1: band pass filtering
firfilter = FIRFilter('hamming',[3.3 , 6.7] , [2.6 , 7.3] , 0.5 , 39.4 , fs);
iirfilter = Cheby2Filter([3.3 , 6.7] , [2.6 , 7.3] , 0.5 , 39.4 , fs);
firfilter.plot(2)
iirfilter.plot(3)

sig1 = firfilter.filter1(obj);
sig2 = iirfilter.filter3(obj);

%% IQ Complex Envelope (phase detection)
midfreq = 0.5*(3.3+6.7);
bw = 5.9-4.1;
sig1iq = sig1.IQ(midfreq , bw);
sig2iq = sig2.IQ(midfreq , bw);

%% Resampling
decimate_factor = 6;
sig1iq.resample(decimate_factor);
sig2iq.resample(decimate_factor);
sig1iq.plot(4)
sig2iq.plot(5)


%% Writing to file
filepathout = fullfile(outputdir , 'signal1_sampled');
sig2iq.writeToFile(filepathout)
