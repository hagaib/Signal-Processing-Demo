clear
close all
disp('Signal 3 Spectral Analysis');

%% Load Signal Data
outputdir = 'output_data';
filename = 'signal3_sampled.mat';
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





%% Noise amplitude and PSD Estimation
% since we filtered and decimated the original signal, the present signal
% has a noise component which is not white anymore, and it's power is
% reduced by a the decimating factor. Taking these into account, it is
% possible to estimate the noise floor using Welch's Method, and average
% the estimated power spectrum in the domain where there is no
% deterministic signal, and the noise wasn't filtered. We can then
% interpolate over the entire spectrum, and multiply by the decimating
% factor to obtain an estimation for the original noise amplitude.

windur =  30;
decimating_factor = 100;
% using a rectangular window to simplify calculations by avoiding noise and
% coherent gain.
[t , f , s] = obj.spectralAnalysis(1, windur, [], 0.5);

disp('Welch method')
M = 5;
N = floor(length(obj.xx)/M);

pwelch = zeros(N , 1);
for i=0:M-1
    x = obj.xx((i*N+1):((i+1)*N));
    sfft= fft(x);
    pwelch = pwelch + sfft .* conj(sfft);
    
end
pwelch = pwelch / M;
figure(2) , plot(fs/N*(0:N-1) , pwelch) , title('Power Spectrum Estimation using Welch''s Method');
noise_amplitude = sqrt(mean(pwelch(25:50))/N*decimating_factor)


%% CW Amplitude Estimation
windur = 3; % ms
[t , f , s] = obj.spectralAnalysis(window_number2, windur, [], 0, 3);
[t , f , s] = obj.spectralAnalysis(window_number1, windur, [], 0, 8);



timeIndex = 10;
x = abs(s(:, timeIndex));
figure(4)
plot(f', x)
title(['DFT amplitude at time ' , num2str(t(timeIndex)) , ' ms, using ' , num2str(windur) , ' ms window']);

[amax,fmax] = extrema(x,1,0);
% The first peak corresponds to the CW spectrum representation
CWf = f(fmax(1));
CWa = amax(1) * 2 / ( floor(windur*fs) * CG1) * 10^(0.1*SL1);
disp(['amplitude estimation is ' num2str(CWa) , ' near frequency ' , num2str(CWf)]) 




%% PMCW Parameter Estimation

bin = 8;
str = ['PMCW component at ' , num2str(f(bin)) , ' kHz bin:'];
disp(str);
figure(5)
x = abs(s(bin,:));
plot(t , x)
title(str)

threshold = 4;
[amax , tmax] = extrema(x,1,1);

tmax = t(tmax(amax>threshold));
amax = amax(amax>threshold);
PRI = mean(tmax(2:end) - tmax(1:end-1))
N = length(amax)
amp_est = mean(amax);
PMCWa = amp_est * 2 / ( floor(windur*fs) * CG1) *10^(0.1*SL1)




%% phase code
% As can be seen from the phase output of the phase detector at 2.5 kHz,
% the phase is not the same between pulses, which means the CW wave is not 
% synchronized with the pulses, A.K.A pulse code 2.
freq = 2.5;
threshold = 0.5;
bandwidth = 0.3;
phasedetector = obj.IQ(freq , bandwidth);
figure(6)
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

windur = 1; % ms
bin =4;
[t , f , s] = obj.spectralAnalysis(window_number1, windur , [] , 0.8);
figure(7)
plot(t , abs(s(bin,:)))
str = ['PMCW component at ' , num2str(f(bin)) , ' kHz bin:'];
title(str)

% using a different window doesn't improve results. Still, tau parameter of
% 2-2.5 can still be estimated.
