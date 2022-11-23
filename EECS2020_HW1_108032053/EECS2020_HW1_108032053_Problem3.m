%% EECS2020 陳凱揚 108032053 Computer HW1 03/12/2021

%% ---------- Codes for Problems 3 ----------
clear all;
close all;
x = double(imread('lena.jpg')); % Lena image, each row is the system input
[M,N]  = size(x); % M: number of rows, N: number of columns
figure
imagesc(x);
colormap(gray); % show gray scale image according to the provided gray-scale colormap
axis image
title("Original Image");
colorbar % show the colormapping

%% ---------- 3(b) ----------
% The original 3(b) system
y = zeros(M, N);
for m = 1:M,
	for n = 3:N-2
		y(m,n) = sum(x(m, n-2:n+2))/5;
	end
end
figure
imagesc(y); 
colormap(gray);
axis image
title("The original 3(b) system");
colorbar

% Re-model the original system into a causal system
y = zeros(M, N);
for m = 1:M,
	for n = 5:N
		y(m,n) = sum(x(m, n-4:n))/5;
	end
end
figure
imagesc(y); 
colormap(gray);
axis image
title("A causal system which be re-modeled from 3(b)")
colorbar

