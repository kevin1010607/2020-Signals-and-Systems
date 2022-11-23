%% EECS2020 陳凱揚 108032053 Computer HW3 05/18/2021

%% ----------
%% ---------- Part 1 Try to be familiar with the properties of the implemented CTFT
%% ----------
% Generate sampled cosine/discrete-time sinusoid
clear all; close all;
F0 = 6;
Fs = 120;
T = 1/Fs;
total_time = 2;
t_axis = (0:T:total_time);
x = cos(2*pi*F0*t_axis);
Npoint = length(x);

% Fourier transform - Analysis
dF = Fs/Npoint;
f_axis = ((1:1:Npoint)-(Npoint+1)/2)*dF;
X = exp(-sqrt(-1)*2*pi*((f_axis.')*t_axis))*(x.')*T;
%% ---------- (a) ----------
% Magnitude spectrum
mag_X = abs(X);   % magnitude
pha_X = angle(X); % phase
% Draw 圖 1-1
figure
plot(f_axis, mag_X,'linewidth',2);
hold
stem(f_axis, mag_X, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('abs(X(F))')
axis([-10 10 -inf inf])
title('Magnitude spectrum')

% Analysis
t = (0:T:6)-2;
sinusoid_full = cos(2*pi*F0*t);
N = length(t)-1;
square_wave = [zeros(1, round(N/3)) ones(1, round(N/3)+1) zeros(1, round(N/3))];
sinusoid_cutoff = sinusoid_full.*square_wave;
% Draw 圖 1-2
set(figure, "position", [200, 50, 1200, 700]);
subplot(3, 1, 1);
plot(t, sinusoid_full);
title("Sinusoid (Full)");
subplot(3, 1, 2);
plot(t, square_wave);
ylabel("Amplitude");
title("Square Wave");
subplot(3, 1, 3);
plot(t, sinusoid_cutoff);
title("Sinusoid (Cutoff)");
xlabel("Time (sec)");

%% ---------- (c) ---------
% Get X(F) by fft()
X_fft = fft(x);
mag_X_fft = abs(X_fft);
f_axis_fft = (0:1:Npoint-1)*dF;
X_c = T*fftshift(X_fft); % shift Fs/2 and multiply T
mag_X_c = abs(X_c);
% Draw 圖 1-3
set(figure, "position", [200, 50, 1200, 700]);
subplot(3, 1, 1);
plot(f_axis, mag_X,'linewidth',2);
hold
stem(f_axis, mag_X, 'r', 'linewidth',1)
title("Magnitude spectrum (X)");
subplot(3, 1, 2);
plot(f_axis_fft, mag_X_fft,'linewidth',2);
hold
stem(f_axis_fft, mag_X_fft, 'r', 'linewidth',1)
ylabel("abs(X(F))");
title("Magnitude spectrum (X\_fft)");
subplot(3, 1, 3);
plot(f_axis, mag_X_c,'linewidth',2);
hold
stem(f_axis, mag_X_c, 'r', 'linewidth',1)
title("Magnitude spectrum (X\_new)");
xlabel("Frequency (kHz)");

%% ---------- (d) ---------
% Zero padding and get new X
x_d = [x zeros(1, 5*(length(t_axis)-1))];
t_axis_d = (0:T:6*total_time);
Npoint_d = length(x_d);
dF_d = Fs/Npoint_d;
f_axis_d = ((1:1:Npoint_d)-(Npoint_d+1)/2)*dF_d;
X_d = T*fftshift(fft(x_d));
mag_X_d = abs(X_d);
% Draw 圖 1-4
figure
plot(t_axis_d, x_d);
xlabel("Time (sec)");
ylabel("Amplitude");
title("Sinusoid (zero padding)");
% Draw 圖 1-5
figure
plot(f_axis_d, mag_X_d,'linewidth',2);
hold
stem(f_axis_d, mag_X_d, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('abs(X(F))')
axis([-10 10 -inf inf])
title('Magnitude spectrum')

%% ---------- (e) ---------
% from -Fs/2 to Fs/2
f_axis_e = ((1:1:Npoint)-(Npoint+1)/2)*dF; 
X_e = zeros(1,length(f_axis_e));
for iFreq = 1:length(f_axis_e)
   for iTime = 1:length(t_axis)
       X_e(iFreq) = X_e(iFreq) + x(iTime)*exp(-sqrt(-1)*2*pi*f_axis_e(iFreq)*t_axis(iTime))*T;
   end
end    
mag_X_e = abs(X_e);
% Draw 圖 1-6
figure
plot(f_axis_e, mag_X_e,'linewidth',2);
hold
stem(f_axis_e, mag_X_e, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('abs(X(F))')
title('Magnitude spectrum [-Fs/2  Fs/2]')

% from -Fs to Fs
f_axis_e = ((1:1:2*Npoint)-(Npoint+1))*dF; 
X_e = zeros(1,length(f_axis_e));
for iFreq = 1:length(f_axis_e)
   for iTime = 1:length(t_axis)
       X_e(iFreq) = X_e(iFreq) + x(iTime)*exp(-sqrt(-1)*2*pi*f_axis_e(iFreq)*t_axis(iTime))*T;
   end
end    
mag_X_e = abs(X_e);
% Draw 圖 1-7
figure
plot(f_axis_e, mag_X_e,'linewidth',2);
hold
stem(f_axis_e, mag_X_e, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('abs(X(F))')
title('Magnitude spectrum [-Fs  Fs]')

% from -2Fs to 2Fs
f_axis_e = ((1:1:4*Npoint)-2*(Npoint+1))*dF; 
X_e = zeros(1,length(f_axis_e));
for iFreq = 1:length(f_axis_e)
   for iTime = 1:length(t_axis)
       X_e(iFreq) = X_e(iFreq) + x(iTime)*exp(-sqrt(-1)*2*pi*f_axis_e(iFreq)*t_axis(iTime))*T;
   end
end    
mag_X_e = abs(X_e);
% Draw 圖 1-8
figure
plot(f_axis_e, mag_X_e,'linewidth',2);
hold
stem(f_axis_e, mag_X_e, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('abs(X(F))')
title('Magnitude spectrum [-2Fs  2Fs]')

%% ---------- (f) ---------
% F0 = 0kHz
F0 = 0;
x = cos(2*pi*F0*t_axis);
f_axis = ((1:1:4*Npoint)-2*(Npoint+1))*dF; 
X = zeros(1,length(f_axis));
for iFreq = 1:length(f_axis)
   for iTime = 1:length(t_axis)
       X(iFreq) = X(iFreq) + x(iTime)*exp(-sqrt(-1)*2*pi*f_axis(iFreq)*t_axis(iTime))*T;
   end
end
mag_X = abs(X);
% Draw 圖 1-9
figure
subplot(2, 1, 1);
stem(0:1:Npoint-1, x,'b', 'linewidth', 2);
xlabel('Time (n)')
ylabel('x[n]');
title('Discrete time sinusoid (F0 = 0kHz)')
subplot(2, 1, 2);
plot(f_axis, mag_X,'linewidth',2);
hold
stem(f_axis, mag_X, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('abs(X(F))')
title('Magnitude spectrum')

% F0 = 6kHz
F0 = 6;
x = cos(2*pi*F0*t_axis);
f_axis = ((1:1:4*Npoint)-2*(Npoint+1))*dF; 
X = zeros(1,length(f_axis));
for iFreq = 1:length(f_axis)
   for iTime = 1:length(t_axis)
       X(iFreq) = X(iFreq) + x(iTime)*exp(-sqrt(-1)*2*pi*f_axis(iFreq)*t_axis(iTime))*T;
   end
end
mag_X = abs(X);
% Draw 圖 1-10
figure
subplot(2, 1, 1);
stem(0:1:Npoint-1, x,'b', 'linewidth', 2);
xlabel('Time (n)')
ylabel('x[n]');
title('Discrete time sinusoid (F0 = 6kHz)')
subplot(2, 1, 2);
plot(f_axis, mag_X,'linewidth',2);
hold
stem(f_axis, mag_X, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('abs(X(F))')
title('Magnitude spectrum')

% F0 = 30kHz
F0 = 30;
x = cos(2*pi*F0*t_axis);
f_axis = ((1:1:4*Npoint)-2*(Npoint+1))*dF; 
X = zeros(1,length(f_axis));
for iFreq = 1:length(f_axis)
   for iTime = 1:length(t_axis)
       X(iFreq) = X(iFreq) + x(iTime)*exp(-sqrt(-1)*2*pi*f_axis(iFreq)*t_axis(iTime))*T;
   end
end
mag_X = abs(X);
% Draw 圖 1-11
figure
subplot(2, 1, 1);
stem(0:1:Npoint-1, x,'b', 'linewidth', 2);
xlabel('Time (n)')
ylabel('x[n]');
title('Discrete time sinusoid (F0 = 30kHz)')
subplot(2, 1, 2);
plot(f_axis, mag_X,'linewidth',2);
hold
stem(f_axis, mag_X, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('abs(X(F))')
title('Magnitude spectrum')

% F0 = 60kHz
F0 = 60;
x = cos(2*pi*F0*t_axis);
f_axis = ((1:1:4*Npoint)-2*(Npoint+1))*dF; 
X = zeros(1,length(f_axis));
for iFreq = 1:length(f_axis)
   for iTime = 1:length(t_axis)
       X(iFreq) = X(iFreq) + x(iTime)*exp(-sqrt(-1)*2*pi*f_axis(iFreq)*t_axis(iTime))*T;
   end
end
mag_X = abs(X);
% Draw 圖 1-12
figure
subplot(2, 1, 1);
stem(0:1:Npoint-1, x,'b', 'linewidth', 2);
xlabel('Time (n)')
ylabel('x[n]');
title('Discrete time sinusoid (F0 = 60kHz)')
subplot(2, 1, 2);
plot(f_axis, mag_X,'linewidth',2);
hold
stem(f_axis, mag_X, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('abs(X(F))')
title('Magnitude spectrum')

% F0 = 90kHz
F0 = 90;
x = cos(2*pi*F0*t_axis);
f_axis = ((1:1:4*Npoint)-2*(Npoint+1))*dF; 
X = zeros(1,length(f_axis));
for iFreq = 1:length(f_axis)
   for iTime = 1:length(t_axis)
       X(iFreq) = X(iFreq) + x(iTime)*exp(-sqrt(-1)*2*pi*f_axis(iFreq)*t_axis(iTime))*T;
   end
end
mag_X = abs(X);
% Draw 圖 1-13
figure
subplot(2, 1, 1);
stem(0:1:Npoint-1, x,'b', 'linewidth', 2);
xlabel('Time (n)')
ylabel('x[n]');
title('Discrete time sinusoid (F0 = 90kHz)')
subplot(2, 1, 2);
plot(f_axis, mag_X,'linewidth',2);
hold
stem(f_axis, mag_X, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('abs(X(F))')
title('Magnitude spectrum')

% F0 = 114kHz
F0 = 114;
x = cos(2*pi*F0*t_axis);
f_axis = ((1:1:4*Npoint)-2*(Npoint+1))*dF; 
X = zeros(1,length(f_axis));
for iFreq = 1:length(f_axis)
   for iTime = 1:length(t_axis)
       X(iFreq) = X(iFreq) + x(iTime)*exp(-sqrt(-1)*2*pi*f_axis(iFreq)*t_axis(iTime))*T;
   end
end
mag_X = abs(X);
% Draw 圖 1-14
figure
subplot(2, 1, 1);
stem(0:1:Npoint-1, x,'b', 'linewidth', 2);
xlabel('Time (n)')
ylabel('x[n]');
title('Discrete time sinusoid (F0 = 114kHz)')
subplot(2, 1, 2);
plot(f_axis, mag_X,'linewidth',2);
hold
stem(f_axis, mag_X, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('abs(X(F))')
title('Magnitude spectrum')

% F0 = 120kHz
F0 = 120;
x = cos(2*pi*F0*t_axis);
f_axis = ((1:1:4*Npoint)-2*(Npoint+1))*dF; 
X = zeros(1,length(f_axis));
for iFreq = 1:length(f_axis)
   for iTime = 1:length(t_axis)
       X(iFreq) = X(iFreq) + x(iTime)*exp(-sqrt(-1)*2*pi*f_axis(iFreq)*t_axis(iTime))*T;
   end
end
mag_X = abs(X);
% Draw 圖 1-15
figure
subplot(2, 1, 1);
stem(0:1:Npoint-1, x,'b', 'linewidth', 2);
xlabel('Time (n)')
ylabel('x[n]');
title('Discrete time sinusoid (F0 = 120kHz)')
subplot(2, 1, 2);
plot(f_axis, mag_X,'linewidth',2);
hold
stem(f_axis, mag_X, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('abs(X(F))')
title('Magnitude spectrum')
