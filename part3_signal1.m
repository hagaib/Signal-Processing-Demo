clear
close all
disp('Signal 1 Spectral Analysis');

%% Load Signal Data
outputdir = 'output_data';
filename = 'signal1_sampled.mat';
filepath = fullfile(fileparts(mfilename('fullpath')) ,outputdir);
filenamepath = fullfile(outputdir , filename);
load(filenamepath)
obj.plot(1)


% Hanning 2 window
window_number1 = 4;
CG1 = 0.5; %coherent gain
SL1 = 1.42; %scallop loss

%Kaiser Bessel 3 window
window_number2 = 26;
CG2 = 0.4; %coherent gain
SL2 = 1.02; %scallop loss

windur = 3; % ms


%% Component Recognition and frequency estimation
%  

%% CW Amplitude Estimation
[t , f , s] = obj.spectralAnalysis(window_number1, windur, [], 0, 2);


timeIndex = 10;
x = abs(s(:, timeIndex));
figure(3)
plot(f', x)
title(['DFT amplitude at time ' , num2str(t(timeIndex)) , ' ms, using ' , num2str(windur) , ' ms window']);

[amax,fmax] = extrema(x,1,1);
% The first peak corresponds to the CW spectrum representation
CWf = 5+ f(fmax(1));
CWa = amax(1) * 2 / ( floor(windur*fs) * CG1);
disp(['amplitude estimation is ' num2str(CWa) , ' near frequency ' , num2str(CWf)]) 
CWf = 5+ fs-f(fmax(2));
CWa = amax(2) * 2 / ( floor(windur*fs) * CG1);
disp(['amplitude estimation is ' num2str(CWa) , ' near frequency ' , num2str(CWf)])


[t , f , s] = obj.spectralAnalysis(window_number2, windur, [], 0, 4);

timeIndex = 10;
x = abs(s(:, timeIndex));
figure(5)
plot(f', x)
title(['DFT amplitude at time ' , num2str(t(timeIndex)) , ' ms, using ' , num2str(windur) , ' ms window']);

[amax,fmax] = extrema(x,1,1);
% The first peak corresponds to the CW spectrum representation
CWf = 5+ f(fmax(1));
CWa = amax(1) * 2 / ( floor(windur*fs) * CG2);
disp(['amplitude estimation is ' num2str(CWa) , ' near frequency ' , num2str(CWf)]) 
CWf = 5+ fs-f(fmax(2));
CWa = amax(2) * 2 / ( floor(windur*fs) * CG2);
disp(['amplitude estimation is ' num2str(CWa) , ' near frequency ' , num2str(CWf)])

% It seems like the Kaiser window's amplitude estimates are less accurate
% than Hanning when white noise is not present. However, the Kaiser window
% detects the weak PMCW much better.


%% PMCW Parameter Estimation
[t , f , s] = obj.spectralAnalysis(window_number1, windur);

bin = 4;
str = ['PMCW component at ' , num2str(f(bin)) , ' kHz bin:'];
disp(str);
figure(6)
x = abs(s(bin,:));
plot(t , x)
title(str)

[amax , tmax] = extrema(x,1,1);
tmax = t(tmax);
PRI = mean(tmax(2:end) - tmax(1:end-1))
N = length(amax)
amp_est = amax(2);
PMCWa = amp_est * 2 / ( floor(windur*fs) * CG1)
% half of the pulses have much smaller amplitude than they should have.
% This is because PRI is 10 but window duration is 3 ms, so the window will
% have good coverage on every second pulse. This can be remedied if window
% overlap is introduced:

[t , f , s] = obj.spectralAnalysis(window_number1, windur , [] , 0.5);
str = ['PMCW component at ' , num2str(f(bin)) , ' kHz bin: , with 50% overlap'];
disp(str);
figure(7)
x = abs(s(bin,:));
plot(t , x)
title(str)

figure
plot(t , angle(s(bin,:)))
% The same analysis is performed with the Kaiser window
[t , f , s] = obj.spectralAnalysis(window_number2, windur , [] , 0.5);

str = ['PMCW component at ' , num2str(f(bin)) , ' kHz bin:'];
disp(str);
figure(8)
x = abs(s(bin,:));
plot(t , x)
title(str)

[amax , tmax] = extrema(x,1,1);
tmax = t(tmax);
PRI = mean(tmax(2:end) - tmax(1:end-1))
N = length(amax)
amp_est = amax(2);
PMCWa = amp_est * 2 / ( floor(windur*fs) * CG2) 

% amplitude estimator is slightly improved


%% phase code
% As can be seen from the phase output of the phase detector at 0.9 kHz,
% the phase remains constant durin pulses, which means there is a single CW
% wave modulated by the pulses, A.K.A pulse code 1.
freq = 0.9;
threshold = 0.1;
bandwidth = 0.3;
phasedetector = obj.IQ(freq , bandwidth);
figure(11)
p1 = subplot(2,1,1);
env = abs(phasedetector.xx);
plot(phasedetector.timevector , env)
str = ['Complex Envelope of singal centered at ' , num2str(freq) , ' kHz'];
title(str);
p2 = subplot(2,1,2);
phase = angle(phasedetector.xx);
phase(env < threshold) = nan;
plot(phasedetector.timevector , unwrap(phase)/(pi))
title('Phase at pulses')
linkaxes([p1,p2] ,'x');



%% tau estimation
% for tau estimation, we need good time resolution. This comes however at
% the expense of frequency resolution. Since we sampled at 3.5 kS/s, for
% 2ms we get an 8 sample window, which means DFT will have a very bad
% frequency resolution.

windur = 2; % ms
bin =3;
[t , f , s] = obj.spectralAnalysis(window_number1, windur , [] , 0.8 , 9);
figure(10)
plot(t , abs(s(bin,:)))
str = ['PMCW component at ' , num2str(f(bin)) , ' kHz bin:'];
title(str)

% using a different window doesn't improve results. Still, tau parameter of
% 2.5-3 can be estimated.

% If we don't want to increase sample rate, we can alternatively look at
% the DFT of the whole signal. The peak height at 0.9 kHz is 34.56. Since
% we found all the other parameters, we can simply calculate:
tau = 34.56/PMCWa*2/PRI/fs
%which is a much better estimation