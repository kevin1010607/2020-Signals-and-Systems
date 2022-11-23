%% EECS2020 陳凱揚 108032053 Computer HW1 03/12/2021

%% ---------- Codes for Problems 2 ----------
clear all;
close all;
f0 = 1024; % frequency in Hz
total_time = 10; % in sec.
fsRatio = 20;
fs = f0*fsRatio; % sampling sampling rate in Hz
T = 1/fs;  % sampling interval in time
A = 1; % magnitude
phi = pi/8; % phase
t = (0:T:total_time);  % time axis
x_DT = A*cos(2*pi*f0*t+phi);
set(figure, "position", [400, 50, 800, 700]);

%% ---------- 2(a) & 2(b) ----------
% fsRatio = 20
fsRatio = 20;
fs = f0*fsRatio;
T = 1/fs;
t = (0:T:total_time);
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 1);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (fsRatio = 20)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]");
sound(x_DT, fs);

% fsRatio = 4
fsRatio = 4;
fs = f0*fsRatio;
T = 1/fs;
t = (0:T:total_time);
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 3);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (fsRatio = 4)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]");
sound(x_DT, fs);

% fsRatio = 2.5
fsRatio = 2.5;
fs = f0*fsRatio;
T = 1/fs;
t = (0:T:total_time);
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 5);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (fsRatio = 2.5)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]");
sound(x_DT, fs);

% fsRatio = 2.2
fsRatio = 2.2;
fs = f0*fsRatio;
T = 1/fs;
t = (0:T:total_time);
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 7);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (fsRatio = 2.2)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]");
sound(x_DT, fs);

% fsRatio = 2
fsRatio = 2;
fs = f0*fsRatio;
T = 1/fs;
t = (0:T:total_time);
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 2);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (fsRatio = 2)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]");
sound(x_DT, fs);

% fsRatio = 1.8
fsRatio = 1.8;
fs = f0*fsRatio;
T = 1/fs;
t = (0:T:total_time);
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 4);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (fsRatio = 1.8)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]");
sound(x_DT, fs);

% fsRatio = 1.4
fsRatio = 1.4;
fs = f0*fsRatio;
T = 1/fs;
t = (0:T:total_time);
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 6);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (fsRatio = 1.4)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]");
sound(x_DT, fs);

% fsRatio = 1.2
fsRatio = 1.2;
fs = f0*fsRatio;
T = 1/fs;
t = (0:T:total_time);
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 8);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (fsRatio = 1.2)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]");
sound(x_DT, fs);

%% ---------- 2(c) ----------
close all;
set(figure, "position", [400, 50, 800, 700]);

% fsRatio = sqrt(2)
fsRatio = sqrt(2);
fs = f0*fsRatio;
T = 1/fs;
t = (0:T:total_time);
x_DT = A*cos(2*pi*f0*t+phi);
subplot(2, 1, 1);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (fsRatio = sqrt(2))');
axis([0, 1/f0*10 -1 1]);
legend("x[n]");

% fsRatio = pi
fsRatio = pi;
fs = f0*fsRatio;
T = 1/fs;
t = (0:T:total_time);
x_DT = A*cos(2*pi*f0*t+phi);
subplot(2, 1, 2);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (fsRatio = \pi)');
axis([0, 1/f0*10 -1 1]);
legend("x[n]");

