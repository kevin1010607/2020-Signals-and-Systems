%%
%   Computer HW1: Listen to sinusoidal signals, perform sampling, and experience system modeling and implementation
%	Sample codes (Matlab script)
%					
%
%
%                                   Edited by Meng-Lin Li, 02/27/2019
%									Modified by Meng-Lin Li, 03/12/2020
%									Modified by Meng-Lin Li, 03/04/2021

%% -----------------------------------------------------
%% ---------- Codes for Problems 1 and 2 ----------
%% ------------------------------------------------------
% ---------- Generate sampled cosine/discrete-time sinusoid ----------
f0 = 1024; % frequency in Hz
total_time = 5; % in sec.

% !!! Sampling in time 
fsRatio = 20;
fs = f0*fsRatio; % sampling sampling rate in Hz
T = 1/fs;  % sampling interval in time
A = 1; % magnitude
phi = pi/8; % phase

t = (0:T:total_time);  % time axis
% x(t) = A*cos(2*pi*f0*t+phi) where A: magnitude (>= 0), f0: 1/(fundamental period), phi: phase
x_CT = A*cos(2*pi*f0*(0:1/(f0*100):total_time)+phi); % x(t), sampling rate is high enough so that x_CT is a good approximation to the CT signal 
x_DT = A*cos(2*pi*f0*t+phi);  % x[n] = x(nT), sampled cosine/discrete time sinusoid
%Npoint = length(x_DT);   % number of points in sampled cosine
sound(x_DT,fs); % Do? reconstruct and play the CT signal by the sound card, type "help sound" (without the double quote) under MATLAB console to see the usage of sound()

figure
% type "help plot" (without the double quote) under MATLAB console to see the usage of plot().
% type "help stem" to see the usage of stem().
%plot((0:1/(f0*100):total_time), x_CT,'-o', 'linewidth', 2);
plot((0:1/(f0*100):total_time), x_CT,'-', 'linewidth', 2); % CT signal
hold on
stem(t, x_DT,'r', 'linewidth', 2); % DT signal
plot(t, x_DT,'r', 'linewidth', 2); % connet the dots, i.e., connect DT x[n]
xlabel('Time (sec.)');
ylabel('x(nT)');
title('Discrete time sinusoid (time domain)');
axis([0 1/f0*5 -A A]); % only observe the signal from time 0 to time 1/f0*3.  once you remove this line, you can see the whole sampled signal (from 0 to 5 sec.)
legend('x(t)', 'x[n]', 'connected x[n]')

% you may try subplot()

% ---------- Problem 1 ----------
% (a) proof x_CT is a periodic signal with a frequency of 1/f0 via MATLAB graphic illustration
% (b) Change A from 0.25 to 1, and tell the changes in the signal you observe and hear
% (c) Change f0 from 256 to 512, to 1024, and tell the changes in the signal you observe and hear
% (d) Change phi from 0 to pi/8, to pi/4, to pi/2, to pi, to 3*pi/2, to 2*pi and tell the changes in the signal you observe and hear, what type of transformation of the independent variable? justify your answer.


% ---------- Problem 2 ----------
% (a) change fsRatio from 20 down to 1.2 (at least try fsRatio = 20, 4, 2.5, 2.2, 2, 1.8, 1.4, and 1.2), and tell any differences among the sounds or any differnces among the DT signals
% (b) following (a), please tell, what feature of the signal, magnitude, frequency, phase
% is changed after the sampling so that you hear the incorrect sound
% (i.e.,the sound is the not same as the sound of x(t))
% (c) Are all the sampled sinusodial signal, i.e., all the DT sinusodial signals, periodic signals?


%% -----------------------------------------------------
%% ---------- Codes for Problem 3 ----------
%% -----------------------------------------------------
x = double(imread('lena.jpg')); % Lena image, each row is the system input
[M,N]  = size(x); % M: number of rows, N: number of columns
figure
% type "help imagesc" (without the double quote) under MATLAB console to see the usage of imagesc().
imagesc(x);
colormap(gray); % show gray scale image according to the provided gray-scale colormap
axis image
colorbar % show the colormapping

%% ---------- 3(a) ----------
y = zeros(M, N-1); % system output
for m = 1:M,
	y(m,:) = x(m, 2:end) - x(m,1:end-1); % What is the system model for the row input? Is this system causal, linear and time invariant? 
end
figure
imagesc(y);
colormap(gray); % show gray scale image according to the provided gray-scale colormap
axis image
colorbar

% See what the system does to the image
figure
imagesc(abs(y)); % supposedly, no negative value for an image; thus take abs() of y.
colormap(gray); % show gray scale image according to the provided gray-scale colormap
axis image
colorbar

%% ---------- 3(b) ----------
y = zeros(M, N); % system output
for m = 1:M,
	for n = 3:N-2
		y(m,n) = sum(x(m, n-2:n+2))/5; % What is the system model for the row input? Is this system causal, linear and time invariant?
	end
end

% See what the system does to the image
figure
imagesc(y); 
colormap(gray); % show gray scale image according to the provided gray-scale colormap
axis image
colorbar

%%
% if the system is not causal, please re-model it into a causal system and implement the causal system you model.
%%

%% ---------- 3(c) ----------
y1 = (x + fliplr(x))/2;
figure
imagesc(y1); 
colormap(gray); % show gray scale image according to the provided gray-scale colormap
axis image
colorbar

y2 = (x - fliplr(x))/2;
figure
imagesc(y2); 
colormap(gray); % show gray scale image according to the provided gray-scale colormap
axis image
colorbar

x_new = y1 + y2; % is x_new equal to x? check it by x-x_new
figure
imagesc(x_new); 
colormap(gray); % show gray scale image according to the provided gray-scale colormap
axis image
colorbar



% ---------- for fun, ----------
% Why is Fs chosen as 44100 Hz?
[x, Fs] = audioread('sister_12sec.wav'); % wavread() for the old version of MATLAB, Assume the sound is recorded at n = 0										
										 % x: the DT audio signal
										 % Fs: sampling rate in Hz
										 										
sound(x,Fs);

y = x(1:2:end); % y[n] = x[2n], time scaling: compression
sound(y,Fs);

y = flipud(x); % y[n] = x[-n], time reversal
sound(y,Fs);
