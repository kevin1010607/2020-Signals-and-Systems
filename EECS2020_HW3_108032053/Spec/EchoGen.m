% 
%   Demof of Echo generation and reverberation
%                                               Edited by Meng-Lin Li,05/02/2019
%                                               Dept. of Electrical Engineering, 
%												National Tsing Hua University, Taiwan

load handel
whos
pause

sound (y, Fs);
pause

% ----- Echo -----
n = 1:length(y);
a = 0.7; % stable

sprintf('Echo/reverberation delay = 100 ms, attenuation factor a = 0.7')
tau = 100e-3;
D = floor(tau*Fs);
ye = filter(1, [1 zeros(1, D-1) -a], y); % comb reverberator
sound (ye, Fs);


%ye2 = filter([-a zeros(1,D-1) 1], [1 zeros(1,D-1) -a], y); % all pass reverberator
%sound(ye2,Fs);


pause
a = 1.1; % unstable
sprintf('Echo/reverberation delay = 200 ms, attenuation factor a = 1.1')
ye = filter(1, [1 zeros(1, D-1) -a], y);
sound(ye, Fs);