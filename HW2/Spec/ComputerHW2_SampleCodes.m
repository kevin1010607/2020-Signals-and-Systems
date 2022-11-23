%%
%   Computer HW2: Convolution sum, echo time estimation, and reverberator implementation
%	sample codes (Matlab script)
%					
%
%
%                                   Edited by Meng-Lin Li, 03/28/2019
%									Modified by Meng-Lin Li, 04/11/2020
%									Modified by Meng-Lin Li, 04/01/2021

clear all; close all;

%% ----------
%% ---------- Part 1 Convolution sum and echo time estimation - noise remover for echo time estimation
%% ----------
% Linear FM/Chirp signal x (if you're interested, you can google "Linear FM signal" or "Chirp signal")
Fs = 100; % Sampling rate in Hz
tau = 10; % Time duration of the linear FM signal in sec
B = 10; % in Hz
beta = B/tau;
A = 1;
t = 0:1/Fs:tau; % time axis
x = A*sin(pi*beta*t.^2); % the linear FM signal

%% --- Hamming windowed linear FM signal, the transmitted sonar signal from the parking sonar
x = hamming(length(x)).'.*x; 

figure
plot(t, x)
xlabel('Time (sec)');
ylabel('Amplitude')
title('Hamming windowed linear FM signal')

%% --- Perfect echo signal without noise, served as an example of how to estimate echo time

y_woNoise = [zeros(1, 1300) x zeros(1,1000)]; % Note that this is not ground truth for Part 1 problems.

% Estimate the echo time from the envelope (what is envelope, see the definition in slide 35, in Topic1_SignalsAndSystems_Part2_HandWriting0317_2021.pdf)
% You can follow the codes here to find out the echo time of noisy y
Envelope = abs(hilbert(y_woNoise)); % Envelope detection, hilbert() converts cosine or sine into its complex exponential counterpart (cos(theta) to exp(jtheta))
[PeakValue, EchoTimeIndex] = max(Envelope); % here we use the peak time as the echo time
EchoTime = EchoTimeIndex*(1/Fs); 

figure
plot( (0:(length(y_woNoise)-1))*(1/Fs), y_woNoise)
xlabel('Time (sec)');
ylabel('Amplitude')
title('Perfect echo signal with its envelope')
hold
plot( (0:(length(y_woNoise)-1))*(1/Fs), Envelope, 'r')

%% ---- Echo signal contaminated by noise
load NoisyEchoSignal % get noisy echo signal y

figure
plot( (0:(length(y)-1))*(1/Fs), y)
xlabel('Time (sec)');
ylabel('Amplitude')
title('Noisy echo signal, Can you find the echo time using your eyes?')

% ---------- (a), (b) ----------
%% ---- Implement your own convolution sum function - MyConv()
% alpha = 0.5;
% x = [1 1 1 1 1];
% h = zeros(1, 7);
% for i = 1:7
%     h(i) = 0.5^(i-1);
% end
% support_x = [0 4];
% support_h = [0 6];
% [y_MyConv, support] = MyConv(x, support_x, h, support_h);
% y_conv = conv(x, h);
% figure
% plot(1:length(y_MyConv), y_MyConv);
% figure
% plot(1:length(y_conv), y_conv);

h = fliplr(x); % impulse response of the noise remover, why do I design the impulse response of the noise remover in this form???
support_h = [0 length(h)-1];
support_y = [0 length(y)-1];
[y_NoiseSuppressed, support_y_NoiseSuppressed] = MyConv(y, support_y, h, support_h); % noise suppressed by the noise remover via your own convolution sum function MyConv
figure
plot(0:length(y_NoiseSuppressed)-1, y_NoiseSuppressed);                                                                    % Function prototype, please see MyConv.m
					
% ---------- (c) ----------
%% ---- Estimate the echo time from the envelope (see the above sample codes)
???

% ---------- (d) ----------
%% what is the minimum achievable error (in seconds) of the estimated echo time??? Justify your answer



%% ----------
%% ---------- Part 2 Reverberator implementation
%% ----------
% ---------- (a) ----------
a = 0.7; % attenuation coef.
D = 5; % digital time delay
UnitImpulse = [1 zeros(1,???)]; % creat unit impulse, starting from n=0;
x = UnitImpulse;
y = filter(1, [1 zeros(1, D-1) -a], x); % help filter, do you know the usage of filter() now? filter([all the bk in order], [all the ak in order], x) where ak and bk are the coef. of the LCCDE in slide 21, in Topic2_LTISystems_Part2 _withNotes.pdf
n = 0:(length(y)-1); % check length of y[n]
figure
stem(n, y); % is it the same as that in your derivation?
xlabel('n')
title('Impulse response')

% ---------- (b) ----------
% Please comment whether this reverberator can be potentially implemented in real time and under what condition of a and D this reverberator is stable. 

% ---------- (c) ----------
% FIR system design and implementation
[x, Fs] = audioread('Halleluyah.wav'); % wavread() for the old version of MATLAB, Assume the sound is recorded from normalized time n = 0
a = 0.7;
tau = 0.1;
D = tau*Fs;
sound(x, Fs);
% FIR approximation
h = [???]; % impulse response of the FIR-approximated reverberator, Note that assume the first sample of h starts from normalized time n = 0;
[y, support_y] = Myconv(???); % what initial condition is assumed when you use your own MyConv()?
                              % If you use MATLAB built in function conv() to verify your own results, help conv(), take care what SHAPE you should use? 'full', 'same', or 'valid', for MATLAB 2018 or up. Give SHAPE a try and see the difference, and Is this an issue we've mentioned in our lecture? (see support and length change after convolution)
sound(y, Fs); % play the created sound with reverberation 
audiowrite('Halleluyah_FIRecho.wav', y, Fs); % wavwrite() for the old version of MATLAB

figure
plot(???)

% ---------- (d) ----------
% IIR system implementation and initial condition
[x, Fs] = audioread('Halleluyah.wav'); % wavread() for the old version of MATLAB
a = 0.7;
D = ?;
sound(x, Fs);
% IIR system implementation (see slide 48, in Topic2_LTISystems_Part2_withNotes.pdf), verify your implementation by MATLAB function filter()
???
y =???
sound(y, Fs); % play the created sound with reverberation 
audiowrite('Halleluyah_IIRecho.wav', y, Fs); % wavwrite() for the old version of MATLAB

figure
plot(???)


% ---------- (e) ----------
% Tuning the attenuation coef. so that you're going to have an unstable system
% What does the output of an unstable reverberator shoulds like?
???

% ---------- (f) ----------
% Tuning the initial condition
???

% ---------- (g) ----------
% Verify the system output is colorized. That is, some frequencies of the sound will be suppressed and some frequencies will be amplified (or emphasized)
% plot the magnitude response
???





