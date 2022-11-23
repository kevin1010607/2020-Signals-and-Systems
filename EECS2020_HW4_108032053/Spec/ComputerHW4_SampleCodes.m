%%
%   Computer HW4:  Discrete-time processing of continuous-time audio signals
%				   (Make your own Karaoke and sound effect)
%
%
%                                   Edited by Meng-Lin Li, 05/28/2020
%									Dept. of Electrical Engineering,
%									National Tsing Hua University, Taiwan
%									Revised by Meng-Lin Li, 05/27/2021
%



%% ---------- FIR filtering by convolution ----------
clear all; close all;

[x, Fs] = audioread('sister_8sec.wav'); % wavread() for the old version of MATLAB, Assume the sound is recorded at n = 0										 
											% x: the DT audio signal
											% Fs: sampling rate in Hz
											


soundsc(x,Fs); % play the music, i.e., reconstruct the music.
             % Note that Values in x are assumed to be in the range -1.0 <= y <= 1.0. Values outside that range will be clipped by sound().
			 % therefore, for this homework, I choose soundsc() to play the music. Please doc soundsc() and sound() to see their difference.

M = length(x); % length of the input signal
t = (0:M-1)*(1/Fs);      % Time axis

%% ---------- Problem 1. Plot the signal in time domain and frequency domain (magnitude spectrum) ----------
%% Can you find the signals of drum beats? What should the signals look like in time domain and in frequency domain?

% for example,
% plot((0:M-1)*(Fs/M), abs(fft(x))); % !!! check the magnitude spectrum, and the frequency ranges from 0 to Fs, dF = Fs/M


%% ---------- Problem 2. Filter the music ----------
% --------- (a) and (b) ----------
Fcut = 4000; % Hz, cut off frequency (you may try 4 kHz cutoff frequency for LPF and 0.5 kHz cutoff for HPF
FilterOrder = 128; % filter order
flags.lowpass = 1; % 1: low pass filter
% low-pass and high-pass filter design using fir1()
% perform filtering using conv()
if flags.lowpass,
	h = fir1(FilterOrder, Fcut/(Fs/2));   % h: impulse response. help or doc fir1(), frequency normalization is done by normalization with Fs/2 instead of Fs in our lectures
else
	h = fir1(FilterOrder, Fcut/(Fs/2), 'high'); % h: impulse response. help fir1(), frequency normalization is done by normalization with Fs/2 instead of Fs in our lectures	
end	

h = h.'; % convert to column vector, because x is a column vector
figure
stem(h); % Plot the impulse response.
figure
freqz(h,1); % Plot the frequency response - log magnitude response and phase response


y = conv(x, h,'same'); % filtering by convolution, with 'same', the length of y will be the same as that of x. Why 'same'? remove the group delay introduced by the FIR LPF filter (think about relationship between group delay and the support change after convolution		
soundsc(y, Fs); %Reconstruction of analog audio signals by the sound card
if flags.lowpass,
    audiowrite(sprintf('Sister_LPF_%dHz.wav',Fcut),y,Fs); % save the processed audio clip
else
    audiowrite(sprintf('Sister_HPF_%dHz.wav',Fcut),y,Fs); % save the processed audio clip
end

% ---------- (c) ----------
% ---------- design a band-stop filter by fir1() ---------
% The goal is to remove human voice while keeping as much of background music as possible.
% You may try cut-off frequencies [100 Hz 4000 Hz] or [80Hz 4000 Hz] for sister_12sec.wav, but the actual
% effectiveness of voice-removal may depend on the pitch range of the
% singer (male/female/children), type of background music, etc.
%
% type 'help fir1' for more details.
%

% ----- design your own band stop filter and perform the filtering -----
% ???
% 


%% FYI
%% Equivalent codes to freqz(h,1)
%% The function of freqz() can be done by the approximated CTFT codes in the computer HW3 or by fft() based on the following codes
Nfft = 512; % zero-padding h up to Nfft points (i.e., sampling the frequency axis with Nfft sample points)
H = fft(h, Nfft); 
dF = Fs/Nfft;
AbsFreqAxis_forFFT = (0:Nfft-1)*dF; %absolute frequency axis for fft(), from 0 to Fs or equivalently from 0 to 2*pi for normalized anaular frequency
NormalizedAngularFreqAxis_forFFT = AbsFreqAxis_forFFT/Fs*2*pi; % normalized angular frequency
figure
subplot(2,1,1)
plot(NormalizedAngularFreqAxis_forFFT, 20*log10(abs(H)));
xlabel('Normalized Angular Frequency')
ylabel('Magnitude (dB)')
axis([0 pi -120 20])
grid on
subplot(2,1,2)
plot(NormalizedAngularFreqAxis_forFFT, unwrap(angle(H))*180/pi); % angle(): find the "principal phase" from -pi to pi, unwrap(): find the continuous phase spectrum
xlabel('Normalized Angular Frequency')
ylabel('Phase (degrees')
axis([0 pi -2500 0])
grid on
%% End of equivalent codes

%% ---------- Problem 3 Re-sampling the ECG and make the ECG audible ----------
load ECG % ECG: ECG signal, ECG_Fs: sampling rate in Hz
figure
plot((0:length(ECG)-1)/ECG_Fs, ECG);
xlabel('Time (in sec.)')
ylabel('Amplitude (a.u.)');
title('ECG signal')

% ---------- (a) ----------
% resampling the ECG signal so that the resampled ECG signal has the same sampling rate as that of the music "sister_8sec.wav"
% Note that after resampling, the resampled ECG has the same data length as x (i.e., the music of sister_8sec.wav)
I = ???; % interpolation ratio
D = ???; % decimation ratio

% your own resampler - composed of "interplator" and "decimator"
???

ECG_resampled = resample(ECG, I, D); % used to verify your own resampler
soundsc(ECG_resampled, Fs); % ECG_Fs*I/D = Fs - the sampling rate of the music. You will find out, you cannot hear the ECG without any further processing.

% ---------- (b) ----------
% figure out why you cannot hear the ECG sound when you simply play the ECG signal with soundsc().
% devise a signal processing procedure (i.e., a DT system) to make the processed ECG audible
% justify your signal processing

???
AudibleECG = ??? % Note that the AudibleECG is with the same sampling as that of the music "sister_8sec.wav" and owns the same data length as the music "sister_8sec.wav"
                 %  

AudibleECG = AudibleECG/max(abs(AudibleECG)); % normalization to the value [-1 1]
soundsc(AudibleECG, Fs);
audiowrite('AudibleECG.wav', AudibleECG, Fs);


% If you'd like to, you can mix the music and the audible ECG.
MixedMusic = AudibleECG + x;
MixedMusic = MixedMusic/max(abs(MixedMusic)); % normalization to the value [-1 1]
soundsc(MxiedMusic, Fs);










