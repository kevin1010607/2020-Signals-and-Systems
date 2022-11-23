%% EECS2020 陳凱揚 108032053 Computer HW4 06/11/2021

%% ----------
%% ---------- Problem 1  Plot the signal in time domain and frequency domain (magnitude spectrum)
%% ----------
clear all; close all;
[x, Fs] = audioread('sister_8sec.wav');
soundsc(x, Fs);
T = 1/Fs;
M = length(x);
t = (0:M-1)*(1/Fs);
dF = Fs/M;
f = ((1:1:M)-M/2)*dF;
X = T*fftshift(fft(x));

% Draw 圖 1-1
set(figure, 'position', [200, 200, 1000, 400]);
plot(t, x);
xlabel('Time (sec)');
ylabel('x(t)');
title('Time domain');
% Draw 圖 1-2
set(figure, 'position', [200, 200, 1000, 400]);
plot(f, abs(X));
xlabel('Frequency (Hz)');
ylabel('abs(X(F))');
title('Frequency domain');

% Listen to a drum beat
soundsc(x(round(2.8*Fs):round(2.84*Fs)), Fs);
% Draw 圖 1-3
set(figure, 'position', [200, 200, 1000, 400]);
plot(t, x);
xlabel('Time (sec)');
ylabel('x(t)');
title('Time domain (Zoom in)');
axis([2.80 2.84 -inf inf]);
% Draw 圖 1-4
set(figure, 'position', [200, 200, 1000, 400]);
plot(f, abs(X));
xlabel('Frequency (Hz)');
ylabel('abs(X(F))');
title('Frequency domain (Zoom in)');
axis([0.5 1.5 0 0.006]);

%% ----------
%% ---------- Problem 2  Filter the music
%% ----------
close all;
%% ---------- (a) ----------
% Low-pass filter
Fcut = 4000;
FilterOrder = [8 16 32 64 128 256];
% 圖 2-1 ~ 圖 2-18
for i = 1:length(FilterOrder)
    h = fir1(FilterOrder(i), Fcut/(Fs/2));
    h = h.';
    y = conv(x, h, 'same');
    Y = T*fftshift(fft(y));
    
    figure
    stem(h);
    xlabel('Sample (n)');
    ylabel('h[n]');
    title(['Impulse response (order = ', num2str(FilterOrder(i)), ')']);
    figure
    freqz(h, 1);
    title(['Frequency response (order = ', num2str(FilterOrder(i)), ')']);
    figure
    plot(f, 20*log10(abs(X)), f, 20*log10(abs(Y)));
    xlabel('Frequency (Hz)');
    ylabel('abs(X(F)) (dB)');
    title(['Spectrum (order = ', num2str(FilterOrder(i)), ')']);
    legend('X(F)', 'Y(F)');
end

%% ---------- (b) ----------
FilterOrder = 256;
% Low-pass filter
Fcut = 4000;
h = fir1(FilterOrder, Fcut/(Fs/2), 'low');
y = conv(x, h, 'same');
soundsc(y, Fs);
% 圖 2-19
figure
freqz(h, 1);
title('Low-pass filter');

% High-pass filter
Fcut = 4000;
h = fir1(FilterOrder, Fcut/(Fs/2), 'high');
y = conv(x, h, 'same');
soundsc(y, Fs);
% 圖 2-20
figure
freqz(h, 1);
title('High-pass filter');

%% ---------- (c) ----------
% Band-stop filter
FilterOrder = 256;
Fcut = [400 9000];
h = fir1(FilterOrder, Fcut/(Fs/2), 'stop');
y = conv(x, h, 'same');
soundsc(y, Fs);
% 圖 2-21
figure
freqz(h, 1);
title('Band-stop filter');

%% ----------
%% ---------- Problem 3  Re-sampling the ECG and make the ECG audible
%% ----------
close all;
load ECG;
%% ---------- (a) ----------
% Resampling
I = 49;
D = 2;
N = length(ECG);
% Sampling rate expander
w = zeros(N+N*(I-1), 1); % Insert I-1 zeros between every two points of ECG signal
idx = 1;
for i = 1:I:length(w)
    w(i) = ECG(idx);
    idx = idx+1;
end
% Low-pass filter
Fcut = min(1/I, 1/D); % normalized frequency, and choose the smaller Fcut to LPF
h = I*fir1(256, Fcut); % Gain = I
v = conv(w, h, 'same');
% Sampling rate compressor
y = zeros(length(v)/D, 1); 
idx = 1;
for i = 1:length(y)
    y(i) = v(idx);
    idx = idx+D;
end

Y = T*fftshift(fft(y));
ECG_resampled = resample(ECG, I, D);
ECG_RESAMPLED = T*fftshift(fft(ECG_resampled));
soundsc(ECG_resampled, Fs);
soundsc(y, Fs);

% Draw 圖 3-3
figure
plot(f, abs(Y));
xlabel('Frequency (Hz)');
ylabel('abs(Y(F))');
title('Frequency domain (Y)');
% Draw 圖 3-4
figure
plot(f, abs(ECG_RESAMPLED));
xlabel('Frequency (Hz)');
ylabel('abs(ECG\_RESAMPLED(F))');
title('Frequency domain (ECG\_RESAMPLED)');
% Draw 圖 3-5
figure
plot(f, abs(Y));
xlabel('Frequency (Hz)');
ylabel('abs(Y(F))');
title('Frequency domain (Y) (Zoom in)');
axis([-200 200 -inf inf]);
% Draw 圖 3-6
figure
plot(f, abs(ECG_RESAMPLED));
xlabel('Frequency (Hz)');
ylabel('abs(ECG\_RESAMPLED(F))');
title('Frequency domain (ECG\_RESAMPLED) (Zoom in)');
axis([-200 200 -inf inf]);

%% ---------- (b) ----------
% Amplitude modulation
Fshift = 1000;
p = cos(2*pi*Fshift*t);
p = p';
r = p.*y;
R = T*fftshift(fft(r));
soundsc(r, Fs);
audiowrite('AudibleECG.wav', r, Fs);

% Draw 圖 3-5
figure
plot(f, abs(R));
xlabel('Frequency (Hz)');
ylabel('abs(R(F))');
title('Frequency domain (R)');
axis([-2000 2000 -inf inf]);

