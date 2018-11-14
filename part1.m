%%
%% Part1 - Signal Generation
%%

clear 
close all

%% Name of input / output directories
inputdir = 'input_data';
outputdir = 'output_data';


%% Data struct to hold signal objects
sigcount = 3;
signals = cell(sigcount,1);

%% Time and sample rate for each signal
time = [100 , 93.5 , 115.5];
fs = [21 , 16 , 1000];

for i=1:sigcount
    filenamein = ['signal' , num2str(i) , '_input.csv'];
    filenameout = ['signal' , num2str(i)];
    filepathin = fullfile(inputdir , filenamein);
    filepathout = fullfile(outputdir , filenameout);
    
    signals{i} = Signal(filepathin);
%     if(i==1), signals{i}.components{3}.I_PHASE_CODE = 2; end
    signals{i}.generate(fs(i) , time(i));
    signals{i}.title = ['Signal ' , num2str(i)];
    signals{i}.writeToFile(filepathout);
end

clear 'filenamein'  'filenameout'  'filepathin'  'filepathout'

clear