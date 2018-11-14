clear
disp('Signal 3 Anti-Aliasing and Resampling');

%% Load Signal Data
outputdir = 'output_data';
filename = 'signal3.mat';
filepath = fullfile(fileparts(mfilename('fullpath')) ,outputdir);
filenamepath = fullfile(outputdir , filename);
load(filenamepath)
obj.plot(1)

%% Signal Anti Aliasing Lowpass Filtering
filter = ButterFilter(3.2 , 3.7 , 0.4 , 34.7 , fs);
filter.plot(2)

sig = filter.filter3(obj);

%% Resampling
decimate_factor = 100;

sig.resample(decimate_factor);
sig.plot(4)


%% Writing to file
filepathout = fullfile(outputdir , 'signal3_sampled');
sig.writeToFile(filepathout)
