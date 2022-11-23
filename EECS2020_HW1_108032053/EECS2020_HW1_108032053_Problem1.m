%% EECS2020 陳凱揚 108032053 Computer HW1 03/12/2021

%% ---------- Codes for Problems 1 ----------
clear all;
f0 = 1024; % frequency in Hz
total_time = 5; % in sec.
fsRatio = 20;
fs = f0*fsRatio; % sampling sampling rate in Hz
T = 1/fs;  % sampling interval in time
A = 1; % magnitude
phi = pi/8; % phase
t = (0:T:total_time);  % time axis

%% ---------- 1(a) ----------
close all;
t = (0:1/(f0*100):total_time);
x_CT = A*cos(2*pi*f0*t+phi);
x_CT_shift = A*cos(2*pi*f0*(t+1/f0)+phi);
difference = x_CT_shift-x_CT;
set(figure, "position", [300, 50, 500, 500]);

% Draw x(t)
subplot(2, 2, 1);
plot(t, x_CT);
xlabel("Time (sec)");
ylabel("y(t)");
title("x(t)");
axis([1/f0 1/f0+1/f0*5 -A A]);
legend("x(t)");

% Draw x(t+1/f0)
subplot(2, 2, 2);
plot(t, x_CT_shift);
xlabel("Time (sec)");
ylabel("y(t)");
title("x(t+1/f0)");
axis([1/f0 1/f0+1/f0*5 -A A]);
legend("x(t+1/f0)");

% Draw x(t) and x(t+1/f0) in the same figure
subplot(2, 2, 3);
plot(t, x_CT, t, x_CT_shift);
xlabel("Time (sec)");
ylabel("y(t)");
title("x(t) & x(t+1/f0)");
axis([1/f0 1/f0+1/f0*5 -A A]);
legend("x(t)", "x(t+1/f0)");

% Draw x(t) - x(t+1/f0)
subplot(2, 2, 4);
plot(t, difference);
xlabel("Time (sec)");
ylabel("y(t)");
title("x(t) - x(t+1/f0)");
axis([1/f0 1/f0+1/f0*5 -A A]);
legend("x(t) - x(t+1/f0)");

%% ---------- 1(b) ----------
close all;
t = (0:T:total_time);
set(figure, "position", [400, 50, 800, 500]);

% A = 1
A = 1;
x_DT = A*cos(2*pi*f0*t+phi);
subplot(2, 1, 1);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
plot(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (A = 1)');
axis([0, 1/f0*3 -1 1]);
legend("x[n]", "connected x[n]");
sound(x_DT, fs);

% A = 0.25
A = 0.25;
x_DT = A*cos(2*pi*f0*t+phi);
subplot(2, 1, 2);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
plot(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (A = 0.25)');
axis([0, 1/f0*3 -1 1]);
legend("x[n]", "connected x[n]");
sound(x_DT, fs);

%% ---------- 1(c) ----------
close all;
A = 1;
set(figure, "position", [400, 50, 800, 700]);

% f0 = 256
f0 = 256;
x_DT = A*cos(2*pi*f0*t+phi);
subplot(3, 1, 1);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
plot(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (f0 = 256)');
axis([0, 1/f0 -1 1]);
legend("x[n]", "connected x[n]");
sound(x_DT, fs);

% f0 = 512
f0 = 512;
x_DT = A*cos(2*pi*f0*t+phi);
subplot(3, 1, 2);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
plot(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (f0 = 512)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]", "connected x[n]");
sound(x_DT, fs);

% f0 = 1024
f0 = 1024;
x_DT = A*cos(2*pi*f0*t+phi);
subplot(3, 1, 3);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
plot(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (f0 = 1024)');
axis([0, 1/f0*4 -1 1]);
legend("x[n]", "connected x[n]");
sound(x_DT, fs);

%% ---------- 1(d) ----------
close all;
f0 = 1024;
set(figure, "position", [400, 50, 800, 700]);

% phi = 0
phi = 0;
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 1);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
plot(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (\phi = 0)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]", "connected x[n]");
sound(x_DT, fs);

% phi = pi/8
phi = pi/8;
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 3);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
plot(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (\phi = \pi/8)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]", "connected x[n]");
sound(x_DT, fs);

% phi = pi/4
phi = pi/4;
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 5);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
plot(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (\phi = \pi/4)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]", "connected x[n]");
sound(x_DT, fs);

% phi = pi/2
phi = pi/2;
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 7);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
plot(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (\phi = \pi/2)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]", "connected x[n]");
sound(x_DT, fs);

% phi = pi
phi = pi;
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 2);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
plot(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (\phi = \pi)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]", "connected x[n]");
sound(x_DT, fs);

% phi = 3*pi/2
phi = 3*pi/2;
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 4);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
plot(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (\phi = 3\pi/2)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]", "connected x[n]");
sound(x_DT, fs);

% phi = 2*pi
phi = 2*pi;
x_DT = A*cos(2*pi*f0*t+phi);
subplot(4, 2, 6);
hold on;
stem(t, x_DT, 'r', 'linewidth', 2);
plot(t, x_DT, 'r', 'linewidth', 2);
xlabel("Time (sec)");
ylabel("x(nT)");
title('Discrete time sinusoid (\phi = 2\pi)');
axis([0, 1/f0*2 -1 1]);
legend("x[n]", "connected x[n]");
sound(x_DT, fs);

