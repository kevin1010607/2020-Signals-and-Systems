%% EECS2020 陳凱揚 108032053 Computer HW3 05/22/2021

%% ----------
%% ---------- Part 2 Apply the implemented CTFT for realistic signal analysis, system design and implementation
%% ----------
clear all; close all;
% Fourier analysis over the given ECG signal
load ECG % ECG: ECG signal, Fs: sampling rate in Hz

figure
plot((0:length(ECG)-1)/Fs, ECG);
xlabel('Time (in sec.)')
ylabel('Amplitude (in mV)');
title('ECG signal')

%% ---------- (a) ----------
T = 1/Fs;
% Single cycle signal
ECG_single = ECG(6610:8064);
Npoint_single = length(ECG_single);
dF_single = Fs/Npoint_single;
f_axis_single = ((1:1:Npoint_single)-(Npoint_single+1)/2)*dF_single;
X_single = T*fftshift(fft(ECG_single));
mag_X_single = abs(X_single);
% Draw 圖 2-1
figure
subplot(2, 1, 1);
plot(((0:1:Npoint_single-1)+6610)/Fs, ECG_single);
xlabel("Time (sec)");
ylabel("x[n]");
title("ECG (single cycle)");
subplot(2, 1, 2);
plot(f_axis_single, mag_X_single);
xlabel("Frequency (Hz)");
ylabel("abx(X(F))");
title("Magnitude spectrum (single cycle)");

% Whole signal
Npoint_whole = length(ECG);
dF_whole = Fs/Npoint_whole;
f_axis_whole = ((1:1:Npoint_whole)-(Npoint_whole+1)/2)*dF_whole;
X_whole = T*fftshift(fft(ECG));
mag_X_whole = abs(X_whole);
% Draw 圖 2-2
figure
subplot(2, 1, 1);
plot((0:1:Npoint_whole-1)/Fs, ECG);
xlabel("Time (sec)");
ylabel("x[n]");
title("ECG (whole)");
subplot(2, 1, 2);
plot(f_axis_whole, mag_X_whole);
xlabel("Frequency (Hz)");
ylabel("abx(X(F))");
title("Magnitude spectrum (whole)");
%% ---------- (b) ----------
% Normalized Amplitude
[max_val_single, max_idx_single] = max(mag_X_single);
[max_val_whole, max_idx_whole] = max(mag_X_whole);
mag_X_single_n = mag_X_single/max_val_single;
mag_X_whole_n = mag_X_whole/max_val_whole;
% Draw 圖 2-3
set(figure, "position", [200, 50, 1200, 700]);
subplot(3, 1, 1);
plot(f_axis_single, mag_X_single);
axis([-60 60 -inf inf]);
ylabel("abs(X(F))");
title("Magnitude spectrum (single cycle)");
subplot(3, 1, 2);
plot(f_axis_whole, mag_X_whole);
axis([-60 60 -inf inf]);
ylabel("abs(X(F))");
title("Magnitude spectrum (whole)");
subplot(3, 1, 3);
plot(f_axis_single, mag_X_single_n, f_axis_whole, mag_X_whole_n);
axis([-60 60 -inf inf]);
xlabel("Frequency (Hz)");
ylabel("Normalized Amplitude");
title("Magnitude spectrum (both)");
legend("single ECG", "whole ECG");

% Analysis
impulse_train = zeros(1, length(ECG));
for i = 900:1500:15000
    impulse_train(i) = 1;
end
% Draw 圖 2-4
set(figure, "position", [200, 50, 1200, 700]);
subplot(3, 1, 1);
plot((0:1:Npoint_single-1)/Fs, ECG_single);
title("ECG (single cycle)");
subplot(3, 1, 2);
plot((0:1:Npoint_whole-1)/Fs, impulse_train);
title("Impulse train");
subplot(3, 1, 3);
plot((0:1:Npoint_whole-1)/Fs, ECG);
xlabel("Time (sec)");
title("ECG (whole)");
%% ---------- (c) ----------
% Draw 圖 2-5
figure
plot(f_axis_whole, mag_X_whole);
axis([-10 10 -inf inf]);
xlabel("Frequency (Hz)");
ylabel("abs(X(F))");
title("Magnitude spectrum (whole)");

%% ---------- (d) ----------
% Design a proper moving average filter
N = 30;
h = [ones(1, N) zeros(1, Npoint_single-N)]*60;
H = T*fftshift(fft(h));
mag_H = abs(H);
% Draw 圖 2-6
figure
plot(f_axis_single, mag_H);
xlabel("Frequency (Hz)");
ylabel("abs(X(F))");
title("Magnitude spectrum (filter)");
% Draw 圖 2-7
figure
plot(f_axis_single, mag_H);
axis([-200 200 -inf inf]);
xlabel("Frequency (Hz)");
ylabel("abs(X(F))");
title("Magnitude spectrum (filter)");

% Time domain
ECG_new = conv(ECG_single, h);
ECG_new = ECG_new(1:1455);
% Draw 圖 2-8
set(figure, "position", [400, 400, 700, 400]);
subplot(2, 1, 1);
plot(((0:1:Npoint_single-1)+6610)/Fs, ECG_single);
xlabel("Time (sec)");
ylabel("x[n]");
title("ECG (before)");
subplot(2, 1, 2);
plot(0:length(ECG_new)-1, ECG_new);
plot(((0:1:Npoint_single-1)+6610)/Fs, ECG_new);
xlabel("Time (sec)");
ylabel("x[n]");
title("ECG (after)");

% Frequency domain
X_new = X_single.*H;
mag_X_new = abs(X_new);
% Draw 圖 2-9
set(figure, "position", [400, 400, 700, 400]);
subplot(1, 2, 1);
plot(f_axis_single, mag_X_single);
axis([-200 200 -inf inf]);
xlabel("Frequency (Hz)");
ylabel("abs(X(F))");
title("Magnitude spectrum (before)");
subplot(1, 2, 2);
plot(f_axis_single, mag_X_new);
axis([-200 200 -inf inf]);
xlabel("Frequency (Hz)");
ylabel("abs(X(F))");
title("Magnitude spectrum (after)");
% Draw 圖 2-10
figure
plot(f_axis_single, mag_X_new, f_axis_single, mag_X_single);
axis([-100 100 -inf inf]);
xlabel("Frequency (Hz)");
ylabel("abs(X(F))");
title("Magnitude spectrum (compare)");
legend("After", "Before");

%% ---------- (e) ----------
% H1(F) = 1/(1-a*exp(-jwD))
a = 0.7;
D = 819;
Npoint = 100000;
dF = 1/Npoint;
f_axis_H1 = ((1:1:Npoint)-(Npoint+1)/2)*dF;
H1 = ones(1, Npoint)./(1-a*exp(-sqrt(-1)*2*pi*f_axis_H1*D));
mag_H1 = abs(H1);
% Draw 圖 2-11
figure
plot(f_axis_H1, mag_H1);
xlabel("Normalized frequency");
ylabel("abs(X(F))");
title("Magnitude spectrum of H1");
% Draw 圖 2-12
figure
plot(f_axis_H1, mag_H1);
axis([0 0.01 0 3.5]);
xlabel("Normalized frequency");
ylabel("abs(X(F))");
title("Magnitude spectrum of H1 (Zoom in)");

% H2(F) = Ye(F)/Y(F)
load handel
n = 1:length(y);
a = 0.7;
T = 1/Fs;
tau = 100e-3;
D = floor(tau*Fs);
Npoint = length(y);
dF = 1/Npoint;
f_axis_H2 = ((1:1:Npoint)-(Npoint+1)/2)*dF;
ye = filter(1, [1 zeros(1, D-1) -a], y);
H2 = fftshift(fft(ye)./fft(y));
mag_H2 = abs(H2);
% Draw 圖 2-13
figure
plot(f_axis_H2, mag_H2);
xlabel("Normalized frequency");
ylabel("abs(X(F))");
title("Magnitude spectrum of H2");
% Draw 圖 2-14
figure
plot(f_axis_H2, mag_H2);
axis([0 0.01 0 3.5]);
xlabel("Normalized frequency");
ylabel("abs(X(F))");
title("Magnitude spectrum of H2 (Zoom in)");

% Compare
% Draw 圖 2-15
figure
plot(f_axis_H1, mag_H1, f_axis_H2, mag_H2);
xlabel("Normalized frequency");
ylabel("abs(X(F))");
title("Magnitude spectrum");
legend("mag_H_1", "mag_H_2");
% Draw 圖 2-16
figure
plot(f_axis_H1, mag_H1, f_axis_H2, mag_H2);
axis([0 0.01 0 3.5]);
xlabel("Normalized frequency");
ylabel("abs(X(F))");
title("Magnitude spectrum (Zoom in)");
legend("mag_H_1", "mag_H_2");

%% ---------- Bonus ----------
% Design a multiple notch filter from the comb reverberator
a = 0.7;
D = 30;
Fs = 1800;
Npoint = length(ECG_new);
dF = 1/Npoint;
f_axis = ((1:1:Npoint)-(Npoint+1)/2)*dF;
H = ones(1, Npoint)./(1-a*exp(-sqrt(-1)*2*pi*f_axis*D));
mag_H = 0.5*(-abs(H)+3.3);
pha_H = angle(H);
% Draw 圖 2-17
figure
plot(f_axis*Fs, mag_H);
axis([-200 200 -inf inf]);
xlabel("Normalized frequency");
ylabel("abs(X(F))");
title("Magnitude spectrum (Zoom in)");
% Draw 圖 2-18
figure
plot(f_axis*Fs, pha_H);
axis([-200 200 -inf inf]);
xlabel("Normalized frequency");
ylabel("angle(X(F))");
title("Phase spectrum (Zoom in)");

% Compare mag_X_comb and mag_X
% Draw 圖 2-19
mag_X_comb = mag_X_single.*mag_H;
figure
plot(f_axis_single, mag_X_comb, f_axis_single, mag_X_single);
axis([-200 200 -inf inf]);
xlabel("Frequency (Hz)");
ylabel("abs(X(F))");
title("Magnitude spectrum (compare Comb and Original)");
legend("mag\_X\_Comb", "mag\_X");

% Compare mag_X_comb and mag_X_MAF (moving average filter)
% Draw 圖 2-20
figure
plot(f_axis_single, mag_X_comb, f_axis_single, mag_X_new);
axis([-200 200 -inf inf]);
xlabel("Frequency (Hz)");
ylabel("abs(X(F))");
title("Magnitude spectrum (compare Comb and MAF)");
legend("mag\_X\_Comb", "mag\_X\_MAF");
















