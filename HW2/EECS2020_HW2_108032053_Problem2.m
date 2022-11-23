a%% EECS2020 陳凱揚 108032053 Computer HW2 04/18/2021

%% ----------
%% ---------- Part 2 Reverberator implementation
%% ----------
clear all; close all;

%% ---------- (a) ----------
% Derive the impulse response of this reverberator by filter()
a = 0.7; % attenuation coef.
D = 5; % digital time delay
M = 10;
UnitImpulse = [1 zeros(1,M*D)]; % creat unit impulse, starting from n=0;
x = UnitImpulse;
y = filter(1, [1 zeros(1, D-1) -a], x);
n = 0:(length(y)-1); % check length of y[n]

% 圖 2-1
% Draw impulse response
figure
stem(n, y);
xlabel('n');
ylabel('Amplitude');
title('Impulse response');
legend('h[n]');

%% ---------- (c) ----------
% FIR system design and implementation
[x, Fs] = audioread('Halleluyah.wav');
a = 0.7;
tau = 0.1; % continuous time delay
D = round(tau*Fs);
sound(x, Fs);

% FIR approximation
% y[n] ~= x[n]+(a)*x[n-D]+(a^2)*x[n-2D]+...+(a^M)*x[n-MD] 
% We ignore next other terms, because a^M < 0.001 --> M = 20
M = 19;
UnitImpulse = [1 zeros(1, M*D)];
h = filter(1, [1 zeros(1, D-1) -a], UnitImpulse);
support_x = [0 length(x)-1];
support_h = [0 length(h)-1];
[y, support_y] = MyConv(x, support_x, h, support_h);
y_conv = conv(x, h);
y_c = y;% for (d)
t = (support_y(1):support_y(2))*(1/Fs);
% Listen and write the signal after change
sound(y, Fs);
audiowrite('Halleluyah_FIRecho.wav', y, Fs);
% Transpose and add zero to x in order to make x and y with same size.
x = [x.' zeros(1, length(y)-length(x))]; 

% 圖 2-2
% Verify by conv()
set(figure, "position", [200, 50, 1200, 700]);
% Draw y(t)
subplot(3, 1, 1);
plot(t, y);
title("y(t)");
legend("y(t)");
% Draw y_conv(t)
subplot(3, 1, 2);
plot(t, y_conv);
ylabel("Amplitude");
title("y_c_o_n_v(t)");
legend("y_c_o_n_v(t)");
% Draw y(t)-y_conv(t)
subplot(3, 1, 3);
plot(t, y-y_conv');
xlabel("Time (sec)");
title("y(t)-y_c_o_n_v(t)");
legend("y(t)-y_c_o_n_v(t)");

% 圖 2-3
% Compare input and output
set(figure, "position", [200, 50, 1200, 700]);
% Draw x(t)
subplot(3, 1, 1);
plot(t, x);
title("x(t)");
legend("x(t)");
% Draw y(t)
subplot(3, 1, 2);
plot(t, y);
ylabel("Amplitude");
title("y(t)");
legend("y(t)");
% Draw y(t)-x(t)
subplot(3, 1, 3);
plot(t, y-x);
xlabel("Time (sec)");
title("y(t)-x(t)");
legend("y(t)-x(t)");

% 圖2-4
% Zoom in y(t)-x(t)
set(figure, "position", [200, 50, 1200, 700]);
subplot(3, 1, 1);
plot(t, y-x);
title("Zoom in");
axis([0.05 0.55 -0.8 0.8]);
legend("y(t)-x(t)");
subplot(3, 1, 2);
plot(t, y-x);
ylabel("y(t)-x(t)");
axis([4.5 5 -0.8 0.8]);
legend("y(t)-x(t)");
subplot(3, 1, 3);
plot(t, y-x);
xlabel("Time (sec)");
axis([8.8 9.3 -0.8 0.8]);
legend("y(t)-x(t)");

%% ---------- (d) ----------
% IIR system implementation and initial condition
[x, Fs] = audioread('Halleluyah.wav');
a = 0.7;
tau = 0.1; % continuous time delay
D = round(tau*Fs);
sound(x, Fs);

% IIR system implementation
y = zeros(1, length(x));
% y[n] = 0, x[n] = 0, when n < 0
for n = 0:D-1
    y(n+1) = 0+x(n+1);
end
for n = D:(length(x)-1)
    y(n+1) = a*y(n+1-D)+x(n+1);
end
y_filter = filter(1, [1 zeros(1, D-1) -a], x);
y_d = y;
t = (0:length(y)-1)*(1/Fs);
% Listen and writethe signal after change
sound(y, Fs);
audiowrite('Halleluyah_IIRecho.wav', y, Fs);

% 圖 2-5
% Verify by filter()
set(figure, "position", [200, 50, 1200, 700]);
% Draw y(t)
subplot(3, 1, 1);
plot(t, y);
title("y(t)");
legend("y(t)");
% Draw y_filter(t)
subplot(3, 1, 2);
plot(t, y_filter);
ylabel("Amplitude");
title("y_f_i_l_t_e_r(t)");
legend("y_f_i_l_t_e_r(t)");
% Draw y(t)-y_filter(t)
subplot(3, 1, 3);
plot(t, y-y_filter');
xlabel("Time (sec)");
title("y(t)-y_f_i_l_t_e_r(t)");
legend("y(t)-y_f_i_l_t_e_r(t)");

% 圖 2-6
% Compare input and output
set(figure, "position", [200, 50, 1200, 700]);
% Draw x(t)
subplot(3, 1, 1);
plot(t, x);
title("x(t)");
legend("x(t)");
% Draw y(t)
subplot(3, 1, 2);
plot(t, y);
ylabel("Amplitude");
title("y(t)");
legend("y(t)");
% Draw y(t)-x(t)
subplot(3, 1, 3);
plot(t, y-x');
xlabel("Time (sec)");
title("y(t)-x(t)");
legend("y(t)-x(t)");

% 圖 2-7
% Compare output of (c) and (d)
set(figure, "position", [200, 50, 1200, 700]);
% Draw y_c(t)
y_c = y_c(1:length(y_d)); % Make y_c and y_d with same size
subplot(3, 1, 1);
plot(t, y_c);
title("y_c(t)");
legend("y_c(t)");
% Draw y_d(t)
subplot(3, 1, 2);
plot(t, y_d);
ylabel("Amplitude");
title("y_d(t)");
legend("y_d(t)");
% Draw y_c(t)-y_d(t)
subplot(3, 1, 3);
plot(t, y_c-y_d);
xlabel("Time (sec)");
title("y_c(t)-y_d(t)");
axis([-inf inf -2 2]);
legend("y_c(t)-y_d(t)");

% 圖 2-8
% Zoom in y_c(t)-y_d(t)
figure
plot(t, y_c-y_d);
xlabel("Time (sec)");
ylabel("Amplitude");
title("y_c(t)-y_d(t)");
legend("y_c(t)-y_d(t)");

%% ---------- (e) ----------
% Tuning the attenuation coef. so that you're going to have an unstable system
% What does the output of an unstable reverberator shoulds like?
[x, Fs] = audioread('Halleluyah.wav');
tau = 0.1;
D = round(tau*Fs);
sound(x, Fs);
t = (0:length(x)-1)*(1/Fs);

% 圖 2-9
% Draw the input ans unstable signal
set(figure, "position", [200, 50, 1200, 700]);
% x(t)
subplot(3, 1, 1);
plot(t, x);
title("x(t)");
legend("x(t)");
% a = 1.05
a = 1.05;
y = filter(1, [1 zeros(1, D-1) -a], x);
sound(y, Fs);
subplot(3, 1, 2);
plot(t, y);
ylabel("Amplitude");
title("a = 1.05");
legend("y(t)");
% a = -1.05
a = -1.05;
y = filter(1, [1 zeros(1, D-1) -a], x);
sound(y, Fs);
subplot(3, 1, 3);
plot(t, y-x);
xlabel("Time (sec)");
title("a = -1.05");
legend("y(t)");

%% ---------- (f) ----------
% Tuning the initial condition
[x, Fs] = audioread('Halleluyah.wav');
a = 0.7;
tau = 0.1; % continuous time delay
D = round(tau*Fs);
sound(x, Fs);
t = (0:length(x)-1)*(1/Fs);

% Calculate and listen
% K = 1
y1 = zeros(1, length(x));
K = 1;
for n = 0:(length(x)-1)
    if n < D
        y1(n+1) = a*K+x(n+1);
    else
        y1(n+1) = a*y1(n+1-D)+x(n+1);
    end
end
sound(y1, Fs);
% K = -1
y2 = zeros(1, length(x));
K = 1;
for n = 0:(length(x)-1)
    if n < D
        y2(n+1) = a*K+x(n+1);
    else
        y2(n+1) = a*y2(n+1-D)+x(n+1);
    end
end
sound(y2, Fs);
% K = a rand number between (-1, 1)
y3 = zeros(1, length(x));
K = -1+2*rand();
for n = 0:(length(x)-1)
    if n < D
        y3(n+1) = a*K+x(n+1);
    else
        y3(n+1) = a*y3(n+1-D)+x(n+1);
    end
end
sound(y3, Fs);

% 圖 2-10
set(figure, "position", [200, 50, 1200, 700]);
% K = 1
subplot(3, 1, 1);
plot(t, y1);
title("K = 1");
legend("y(t)");
% K = -1
subplot(3, 1, 2);
plot(t, y2);
ylabel("Amplitude");
title("K = -1");
legend("y(t)");
% K = a rand number between (-1, 1)
subplot(3, 1, 3);
plot(t, y3);
xlabel("Time (sec)");
title("K = a rand number between(-1, 1)");
legend("y(t)");

% 圖 2-11
% Zoom in
set(figure, "position", [200, 50, 1200, 700]);
% K = 1
subplot(3, 1, 1);
plot(t, y1);
title("K = 1");
axis([0 2 -2 2]);
legend("y(t)");
% K = -1
subplot(3, 1, 2);
plot(t, y2);
ylabel("Amplitude");
title("K = -1");
axis([0 2 -2 2]);
legend("y(t)");
% K = a rand number between (-1, 1)
subplot(3, 1, 3);
plot(t, y3);
xlabel("Time (sec)");
title("K = a rand number between(-1, 1)");
axis([0 2 -2 2]);
legend("y(t)");

%% ---------- (g) ----------
% Verify the system output is colorized. That is, some frequencies of the sound will be suppressed and some frequencies will be amplified (or emphasized)
% plot the magnitude response
A = 1;
phi = pi/8;
a = 0.7;
tau = 0.1;

% Calculate
% Fs = 50 → D = 5
Fs = 50;
D = round(tau*Fs); 
f0 = 0:(Fs/2)/10000:(Fs/2);
t = 0:1/Fs:1;
AmplitudeRatio1 = zeros(1, 10001);
idx = 1;
for i = f0
    x = A*cos(2*pi*i*t+phi);
    y = filter(1, [1 zeros(1, D-1) -a], x);
    AmplitudeRatio1(idx) = sum(abs(y))/sum(abs(x));
    idx = idx+1;
end
% Fs = 100 → D = 10
Fs = 100;
D = round(tau*Fs);
f0 = 0:(Fs/2)/10000:(Fs/2);
t = 0:1/Fs:1;
AmplitudeRatio2 = zeros(1, 10001);
idx = 1;
for i = f0
    x = A*cos(2*pi*i*t+phi);
    y = filter(1, [1 zeros(1, D-1) -a], x);
    AmplitudeRatio2(idx) = sum(abs(y))/sum(abs(x));
    idx = idx+1;
end
% Fs = 200 → D = 20;
Fs = 200;
D = round(tau*Fs);
f0 = 0:(Fs/2)/10000:(Fs/2);
t = 0:1/Fs:1;
AmplitudeRatio3 = zeros(1, 10001);
idx = 1;
for i = f0
    x = A*cos(2*pi*i*t+phi);
    y = filter(1, [1 zeros(1, D-1) -a], x);
    AmplitudeRatio3(idx) = sum(abs(y))/sum(abs(x));
    idx = idx+1;
end

% 圖 2-12
% Draw
set(figure, "position", [200, 50, 1200, 700]);
% Fs = 50 → D = 5
subplot(3, 1, 1);
plot(f0/Fs, AmplitudeRatio1);
set(gca, 'xtick', [0 0.125 0.25 0.375 0.5]);
set(gca, 'xticklabel', {'0', 'Fs/8', 'Fs/4', '3Fs/8', 'Fs/2'});
title("Fs = 50 → D = 5");
legend("abs(A'/A)");
% Fs = 100 → D = 10
subplot(3, 1, 2);
plot(f0/Fs, AmplitudeRatio2);
set(gca, 'xtick', [0 0.125 0.25 0.375 0.5]);
set(gca, 'xticklabel', {'0', 'Fs/8', 'Fs/4', '3Fs/8', 'Fs/2'});
ylabel("abs(A'/A)");
title("Fs = 100 → D = 10");
legend("abs(A'/A)");
% Fs = 200 → D = 20
subplot(3, 1, 3);
plot(f0/Fs, AmplitudeRatio3);
set(gca, 'xtick', [0 0.125 0.25 0.375 0.5]);
set(gca, 'xticklabel', {'0', 'Fs/8', 'Fs/4', '3Fs/8', 'Fs/2'});
xlabel("frequency");
title("Fs = 200 → D = 20");
legend("abs(A'/A)");
