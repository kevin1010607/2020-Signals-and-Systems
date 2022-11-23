%% EECS2020 陳凱揚 108032053 Computer HW2 04/18/2021

%% ----------
%% ---------- Part 1 Convolution sum and echo time estimation - noise remover for echo time estimation
%% ----------
% Linear FM/Chirp signal x
clear all; close all;
Fs = 100; % Sampling rate in Hz
tau = 10; % Time duration of the linear FM signal in sec
B = 10; % in Hz
beta = B/tau;
A = 1;
t = 0:1/Fs:tau; % time axis
x = A*sin(pi*beta*t.^2); % the linear FM signal
x = hamming(length(x)).'.*x;
figure
plot(t, x)
xlabel('Time (sec)');
ylabel('Amplitude')
title('Hamming windowed linear FM signal')

% Echo signal contaminated by noise
load NoisyEchoSignal % get noisy echo signal y
figure
plot((0:(length(y)-1))*(1/Fs), y)
xlabel('Time (sec)');
ylabel('Amplitude')
title('Noisy echo signal, Can you find the echo time using your eyes?')

%% ---------- (b) ----------
% ---- Implement your own convolution sum function - MyConv()

% Example 2.4 (slide 50, Topic2_LTISystems_Part1_HandWriting0331_2021.pdf)
alpha = 2;
u = ones(1, 5);
h = ones(1, 7);
for i = 2:7
    h(i) = alpha*h(i-1);
end
support_u = [0 4];
support_h = [0 6];
[v_MyConv, support_v] = MyConv(u, support_u, h, support_h);
v_conv = conv(u, h);

% 圖 1-4
set(figure, "position", [200, 50, 1200, 700]);
% Draw v_MyConv[n]
subplot(3, 1, 1);
plot(support_v(1):support_v(2), v_MyConv);
title("v_M_y_C_o_n_v[n]");
legend("v_M_y_C_o_n_v[n]");
% Draw v_conv[n]
subplot(3, 1, 2);
plot(0:(length(v_conv)-1), v_conv);
ylabel("Amplitude");
title("v_c_o_n_v[n]");
legend("v_c_o_n_v[n]");
% Draw v_MyConv[n]-v_conv[n]
subplot(3, 1, 3);
plot(support_v(1):support_v(2), v_MyConv-v_conv);
xlabel("Sample");
title("v_M_y_C_o_n_v[n] - v_c_o_n_v[n]");
legend("v_M_y_C_o_n_v[n] - v_c_o_n_v[n]");

%% ---------- (c) ----------
% impulse response of the noise remover, why do I design the impulse response of the noise remover in this form???
h = fliplr(x);
support_h = [0 length(h)-1];
support_y = [0 length(y)-1];
% noise suppressed by the noise remover via your own convolution sum function MyConv
[y_NoiseSuppressed, support_y_NoiseSuppressed] = MyConv(y, support_y, h, support_h);
t = (support_y_NoiseSuppressed(1):support_y_NoiseSuppressed(2))*(1/Fs);
% Envelope detection, hilbert() converts cosine or sine into its complex exponential counterpart (cos(theta) to exp(jtheta))
Envelope = abs(hilbert(y_NoiseSuppressed));

% 圖 1-5
% Draw y_NoiseSuppressed and Envelope
set(figure, "position", [200, 50, 1200, 700]);
plot(t, y_NoiseSuppressed); 
xlabel("Time(sec)");
ylabel("Amplitude");
title("y_N_o_i_s_e_S_u_p_p_r_e_s_s_e_d");
hold
plot(t, Envelope, 'r');
legend("y_N_o_i_s_e_S_u_p_p_r_e_s_s_e_d[n]", "Envelope");

% 圖 1-6
% Zoom in
set(figure, "position", [200, 50, 1200, 700]);
plot(t, y_NoiseSuppressed); 
xlabel("Time(sec)");
ylabel("Amplitude");
title("y_N_o_i_s_e_S_u_p_p_r_e_s_s_e_d(Zoom in)");
axis([23 27 -150 200]);
hold
plot(t, Envelope, 'r');
legend("y_N_o_i_s_e_S_u_p_p_r_e_s_s_e_d[n]", "Envelope");


% Estimate the echo time from the envelope
[PeakValue, PeakIndex] = max(Envelope);
DelayIndex = (length(h)-1)/2;
EchoTimeIndex = PeakIndex-DelayIndex;
EchoTime = EchoTimeIndex*(1/Fs);

