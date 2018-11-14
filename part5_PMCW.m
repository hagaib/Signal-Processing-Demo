SNRin = [0 , -5 , -10 , -20 , -30 , -40];
%% Change this parameter to set min SNR value for simulation
Len = length(SNRin)-1;

A = 10.^(SNRin/20);
w = 1.0; %kHz
k = [1000 , 1000.5]; % bin number (10.5 for in between bins)
PRI = 100;
tau = 10;
pulsecount = 10;

%% Data for Rectanglular and Kaiser Bessel 2 windows:
WCPL = [3.92 , 3.20];
PL = [1 , 1.5];


%% Desired SNR output at kth DFT bin:
SNRout = 14;



PLdb = 10*log10(PL);

N_rect_on = floor(4*(PRI/tau)^2*10.^(0.1*(SNRout-SNRin+PLdb(1))));
N_rect_mid = floor(4*(PRI/tau)^2*10.^(0.1*(SNRout-SNRin+WCPL(1))));
N_kaiser_on = floor(4*(PRI/tau)^2*10.^(0.1*(SNRout-SNRin+PLdb(2))));
N_kaiser_mid = floor(4*(PRI/tau)^2*10.^(0.1*(SNRout-SNRin+WCPL(2))));


% For Convenience use even number for N
N_rect_on = N_rect_on + mod(N_rect_on,2);
N_rect_mid = N_rect_mid + mod(N_rect_mid,2);
N_kaiser_on = N_kaiser_on + mod(N_kaiser_on,2);
N_kaiser_mid = N_kaiser_mid + mod(N_kaiser_mid,2);


varnames = {'N_rect_on' , 'N_rect_mid' , 'N_kaiser_on' , 'N_kaiser_mid'};

%% CW DATA on bin (rect win)
% CW DFT peak: A/2*N
% noise mean: sqrt(N)
% output SNR: 20*log10(A/2*sqrt(N))

%% CW Data on bin (any window)
% output SNR: 10*log10((A/2)^2*N * PG)

%% CW Data off bin (any window)
% output SNR: 10*log10(N*(A/2)^2) - WCPL


wingen = WindowGenerator();


%% Rect window fixed on a bin k
N = N_rect_on;
fs = N / k(1) * w;

disp('Rectangular window on the bin');
for i=1:Len
    signal = Signal([1 , 1]);
    signal.addComp([3 , A(i) , w , tau , PRI , pulsecount , 1]);
    signal.generate(fs(i) , N(i)/fs(i));
    win = wingen.generate(1, N(i)+1 , 'dfteven');
    plot((0:N(i)-1)*fs(i)/N(i), abs(fft(win.*signal.xx)));
    title(['CW signal in noise ,  SNR = ', num2str(SNRin(i)) ,  ' , N = ' , num2str(N(i))]);
    snr_bin = 10*log10((A(i)/2*tau/PRI)^2*N(i) / PL(1));
    disp(['Amp SNR:' , num2str(SNRin(i)) , ' , calculated SNR in bin:' , num2str(snr_bin) , ' , N = ' , num2str(N(i))]);
end

%% Rect window fixed midway between bins 
N = N_rect_mid;
fs = N / k(2) * w;

disp('Rectangular window midway between bins');
for i=1:Len
    signal = Signal([1 , 1]);
    signal.addComp([3 , A(i) , w , tau , PRI , pulsecount , 1]);
    signal.generate(fs(i) , N(i)/fs(i));
    win = wingen.generate(1, N(i)+1 , 'dfteven');
    plot((0:N(i)-1)*fs(i)/N(i), abs(fft(win.*signal.xx)));
    title(['CW signal in noise ,  SNR = ', num2str(SNRin(i)) ,  ' , N = ' , num2str(N(i))]);
    snr_bin = 10*log10((A(i)/2*tau/PRI)^2*N(i)) - WCPL(1);
    disp(['Amp SNR:' , num2str(SNRin(i)) , ' , calculated SNR in bin:' , num2str(snr_bin) , ' , N = ' , num2str(N(i))]);
end


%% Kaiser window fixed on a bin k
N = N_kaiser_on;
fs = N / k(1) * w;

disp('Kaiser window on the bin');
for i=1:Len
    signal = Signal([1 , 1]);
    signal.addComp([3 , A(i) , w , tau , PRI , pulsecount , 1]);
    signal.generate(fs(i) , N(i)/fs(i));
    win = wingen.generate(24, N(i)+1 , 'dfteven');
    plot((0:N(i)-1)*fs(i)/N(i), abs(fft(win.*signal.xx)));
    title(['CW signal in noise ,  SNR = ', num2str(SNRin(i)) ,  ' , N = ' , num2str(N(i))]);
    snr_bin = 10*log10((A(i)/2*tau/PRI)^2*N(i) / PL(2));
    disp(['Amp SNR:' , num2str(SNRin(i)) , ' , calculated SNR in bin:' , num2str(snr_bin) , ' , N = ' , num2str(N(i))]);
end

%% Rect window fixed midway between bins 
N = N_kaiser_mid;
fs = N / k(2) * w;

disp('Kaiser window midway between bins');
for i=1:Len
    signal = Signal([1 , 1]);
    signal.addComp([3 , A(i) , w , tau , PRI , pulsecount , 1]);
    signal.generate(fs(i) , N(i)/fs(i));
    win = wingen.generate(24, N(i)+1 , 'dfteven');
    plot((0:N(i)-1)*fs(i)/N(i), abs(fft(win.*signal.xx)));
    title(['CW signal in noise ,  SNR = ', num2str(SNRin(i)) ,  ' , N = ' , num2str(N(i))]);
    snr_bin = 10*log10((A(i)/2*tau/PRI)^2*N(i)) - WCPL(2);
    disp(['Amp SNR:' , num2str(SNRin(i)) , ' , calculated SNR in bin:' , num2str(snr_bin) , ' , N = ' , num2str(N(i))]);
end


