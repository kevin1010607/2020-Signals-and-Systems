%
%   Computer HW3:  Experiencing your first Fourier analysis and filter design using a computer
%	sample codes (Matlab script - Example of CTFT and ICTFT implementation)
%					
%                                   Edited by Meng-Lin Li, 05/02/2019
%									Revised by Meng-Lin Li, 05/07/2020
%									Revised by Meng-Lin Li, 05/06/2021
%									Dept. of Electrical Engineering,
%									National Tsing Hua University, Taiwan
%

%% ---------- Part 1 ----------
%% ---------- Generate sampled cosine/discrete-time sinusoid ----------
F0 = 6; % in kHz
Fs = 120; % sampling rate/sampling frequency, in kHz or ksamples/sec
T = 1/Fs;  % time resolution, i.e., sampling interval in time domain, in ms
total_time = 2; % in ms

% !!! Sampling in time 
t_axis = (0:T:total_time);  % time axis
x = cos(2*pi*F0*t_axis);  % sampled cosine/discrete time sinusoid ,time domain
Npoint = length(x);   % number of points in sampled cosine

figure
plot(t_axis, x,'linewidth',2);
hold
stem(t_axis, x,'r','linewidth',2);
xlabel('Time (ms)');
ylabel('x(nT)');
title('Discrete time sinusoid (time domain)');

figure
stem(0:1:Npoint-1, x,'b', 'linewidth', 2);
xlabel('Time (n)')
ylabel('x[n]');
title('Discrete time sinusoid (time domain)')

%% ---------- Fourier transform - Analysis ----------
% !!! Sampling in frequency
dF = Fs/Npoint; % frequency resolution, i.e., sampling interval in frequency domain
%f_axis = (0:1:(Npoint-1))*dF;   % frequency axis (from 0 to Fs or equivalently from 0 to 2*pi for normalized angular frequency)
%f_axis = ((1:1:Npoint)-(Npoint+1)/2)*dF; % frequency axis in kHz (from -Fs/2 to Fs/2 or equivalently from -pi to +pi for normalized angular frequency)
f_axis = ((1:1:Npoint)-(Npoint+1)/2)*dF; % frequency axis in kHz (from -Fs/2 to Fs/2 or equivalently from -pi to +pi for normalized angular frequency)
X = zeros(1,length(f_axis)); % spectrum

% implementatoin of X(Fk) = summation x(nT)*exp(-j*2*pi*Fk*(nT))*T 
for iFreq = 1:length(f_axis),
    iFreq
   for iTime = 1:length(t_axis),
       X(iFreq) = X(iFreq) + x(iTime)*exp(-sqrt(-1)*2*pi*f_axis(iFreq)*t_axis(iTime))*T; % what if t_axis is not starting from time = 0?
   end
end    

%% Part 1 (c)
% !!! You can compare the result with that from MATLAB fft()
%X = fft(x); % spectrum of sampled cosine, frequency domain, complex; Any difference from the spectrum X obtained by our own CTFT codes?
%f_axis_forFFT = ???; %frequency axis for fft(x), from 0 to Fs or equivalently from 0 to 2*pi for normalized angular frequency
                      % with the help of MATLAB fftshift(), frequency axis will become from -Fs/2 to Fs/2
 
mag_X = abs(X);   % magnitude
pha_X = angle(X); % phase

figure
plot(f_axis, mag_X,'linewidth',2);
hold
stem(f_axis, mag_X, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('abs(X(F))')
title('Magnitude spectrum')

figure
plot(f_axis, pha_X,'linewidth',2);
hold
stem(f_axis, pha_X, 'r', 'linewidth',1)
xlabel('Frequency (kHz)');
ylabel('phase(X(F))')
title('Phase spectrum')

% figure
% subplot(2,1,1)
% plot(t_axis, x,'linewidth',2);
% hold
% stem(t_axis, x,'r','linewidth',2);
% axis([0.2 0.8 -1 1]);   % to zoom in
% set(gca,'fontsize',14);
% set(gca,'linewidth',2);
% set(gca,'box','off');
% xlabel('Time (ms)');
% title('Sampled cosine (time domain): x(nT)');
% legend('Original cosine','Sampled cosine','0');	
% legend('boxoff')
% 
% subplot(2,1,2)
% plot(f_axis, mag_X,'linewidth',2);
% set(gca,'fontsize',14);
% set(gca,'linewidth',2);
% set(gca,'box','off');
% xlabel('Frequency (kHz)');
% title('Magnitude spectrum (frequency domain)')
% shg


%% ---------- Inverse Fourier transform - Synthesis (For your information if you're interested) ----------
xx = zeros(1, length(t_axis));
for iTime = 1:length(t_axis),
	iTime
	for iFreq = 1:length(f_axis),
		xx(iTime) = xx(iTime) + X(iFreq)*exp(sqrt(-1)*2*pi*f_axis(iFreq)*t_axis(iTime))*dF;
	end
end
% --- Mean Squared Errors ---
MSE = sum(abs(xx-x).^2)

figure
subplot(2,1,1)
plot(t_axis, x, 'r', 'linewidth',2);
hold
plot(t_axis, real(xx), 'k-.', 'linewidth',2); % imag(xx) should be close to 0
axis([0.2 0.8 -1 1]);   % to zoom in
set(gca,'fontsize',14);
set(gca,'linewidth',2);
set(gca,'box','off');
xlabel('Time (ms)');
legend('Original cosine','Synthesized cosine','0');	% new !!!!!
legend('boxoff')

subplot(2,1,2)
plot(t_axis, abs(xx-x), 'linewidth',2);
xlabel('Time (ms)');
title('Trashogram');
set(gca,'fontsize',14);
set(gca,'linewidth',2);

%% ---------- Part 2 ---------- 
%% (a) 
% Fourier analysis over the given ECG signal
load ECG % PPG: PPG signal, Fs: sampling rate in Hz

figure
plot((0:length(ECG)-1)/Fs, ECG);
xlabel('Time (in sec.)')
ylabel('Amplitude (in mV)');
title('ECG signal')

%% (b) 
% Tell the difference between Fourier analysis over (1) one single-cycle ECG wavelet (i.e., one heart-beat cycle) and (2) the whole given ECG signal

%% (c)
% From the Fourier analysis of the given ECG signal, can you figure out the heart rate?

%% (d)
% power line noise-reduction filter design

%% (e)
% magnitude response of the comb reverberator 









