clear
close all
disp('Signal 2 Spectral Analysis');

%% Load Signal Data
outputdir = 'output_data';
filename = 'signal2_sampled.mat';
filepath = fullfile(fileparts(mfilename('fullpath')) ,outputdir);
filenamepath = fullfile(outputdir , filename);
load(filenamepath)
obj.plot(1)


%Hamming window
window_number = 7;
CG = 0.54; % coherent gain
SL = 1.78; % scallop loss

%% Component 1 - CW wave
% For the CW component, it is better to take a longer window because of the
% trade off between time and frequency resolution introduced by the
% uncertainty principle
disp('CW component')
windur = 10; % ms
[t , f , s] = obj.spectralAnalysis(window_number, windur, [], 0, 2);
% A thin horizontal line at the 0.4 kHz bin appears.

%Let's check for amplitude:
timeIndex = floor(size(s,2)/2);
x = abs(s(:, timeIndex));
figure(3)
plot(f', x)
title(['DFT amplitude at time ' , num2str(t(timeIndex)) , ' ms, using ' , num2str(windur) , ' ms window']);


[amax,fmax] = extrema(x,1,0);
% The first peak corresponds to the CW spectrum representation
disp('frequency estimation')
CWf = f(fmax(1))
CWa = amax(1);
disp('amplitude estimation')
CWa = amax(1) * 2 / ( floor(windur*fs) * CG)
disp('amplitude estimation with scallop loss')
CWa = amax(1) * 2 / ( floor(windur*fs) * CG)*10^(SL/20)



%% Components 2 and 3 - Pulse Modulated CW wave
disp('Pulse Modulated CW Component 1');
% better time resolution is required. Lets take a shorter window:
windur = 2; % ms
[t , f , s] = obj.spectralAnalysis(window_number, windur,[],0, 5);

% 2 pulse trains appear: strong pulses around ~ 1.4kHz and weaker pulses
% around 2.5 kHz.
%Checking for better frequency estimation by looking at figure 3, we see
%that the strong pulse corresponds to the peak at 1.2 kHz, and the weak one
%at 2.2 kHz. (which is close)


%% amplitude and PRI and pulse count estimation:
bin = 3;
str = ['PMCW component at ' , num2str(f(bin)) , ' kHz bin:'];
disp(str);
figure(6)
x = abs(s(bin,:));
plot(t , x)
title(str)
% 8 pulses are clearly visible
[amax , tmax] = extrema(x,1,0);
tmax = t(tmax);
threshold = length(obj.xx)/55;
realpeaks = find(amax>threshold);
amax = amax(realpeaks);
tmax = tmax(realpeaks);
PRI = mean(tmax(2:end) - tmax(1:end-1))
N = length(amax)
amp_est = mean(amax);
PMCWa = amp_est * 2 / ( floor(windur*fs) * CG)*10^(SL/20)


%% amplitude and PRI and pulse count estimation:
bin = 6;
str = ['PMCW component at ' , num2str(f(bin)) , ' kHz bin:'];
disp(str);
threshold = length(obj.xx)/55;
figure(7)
x = abs(s(bin,:));
plot(t , x)
title(str);
% 11 pulses are clearly visible
[amax , tmax] = extrema(x,1,1);
tmax = t(tmax);
realpeaks = find(amax>threshold);
amax = amax(realpeaks);
tmax = tmax(realpeaks);
PRI = mean(tmax(2:end) - tmax(1:end-1))
N = length(amax)
amp_est = mean(amax);
PMCWa = amp_est * 2 / ( floor(windur*fs) * CG)*10^(SL/20)


%% tau estimation
% for tau estimation, we need good time resolution. This comes however at
% the expanse of frequency resolution. Since we sampled at 8 kS/s, for
% 0.5ms we get 4 sample window, which means DFT bins will be 2kHz wide.
% Because of spectral leakage, we will not have enough frequency resolution to view each component
% separately.
windur = 0.5; % ms
[t , f , s] = obj.spectralAnalysis(window_number, windur , [] , 0);
figure(8)
plot(t , abs(s(end,:)))
title(['2 PMCW components on the same frequency bin, using a ', num2str(windur), 'ms Hamming window']);


% Can this be improved using an exact Blackman window?
[t , f , s] = obj.spectralAnalysis(28, windur , [] , 0 , 9);
figure(10)
plot(t , abs(s(end,:)))
title(['2 PMCW components on the same frequency bin, using a ', num2str(windur), 'ms Exact-Blackman window']);
% Not really


% By inspecting only the separated peaks, measuring peak length, we arrive
% at approximately 2 ms (5 sample points) for each. adding twice the transition
% time of the window from 0 to 1 which is half te length, we arrive at 2 +
% 2*(0.25) = 2.5 ms, which is the exact value for tau


%% phase code
% We find the complex envelope at both 1.2 kHz and 2.3 frquencies, and
% plotting both envelope and phase, we see that both components exhibit
% constant phase during pulses. However, the 1.2 kHz component' phase changes
% every pulse, while the 2.3 kHz component's phase is approximately the same 
% in every pulse. From this we can conclude that the
% 2.3 khz component's CW is not dependent on the pulses, so it's not synchronized 
% with the pulses (phase code 1), while the 1.2 kHz component's CW is synchronized 
%with pulses (phase code 2).

freq = 1.2;
threshold = 1.5;
bandwidth = 0.4;
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

freq = 2.3;
threshold = 1.5;
bandwidth = 0.25;
phasedetector = obj.IQ(freq , bandwidth);
figure(12)
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
% figure(3)
% plot(f,abs(s(:,10)))
% mean(abs(s(2,:)))
% floor(2.9*fs)*2.3/2*0.54*10^(-1.78/20)
