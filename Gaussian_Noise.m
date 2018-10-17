%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% Gaussian_Noise.m: Generate simulated Gaussian noise

%% Produce Gaussian noise
D = 500;                                % Data set size
I = randn(1,D);                         % In-phase component
Q = randn(1,D);                         % Quadrature component
GaussianNoise = (I + 1i.*Q)/sqrt(2);	% Produce 1xD matrix of Gaussian noise

% Plot data
figure;
n = 0:(D - 1);
plot(n,GaussianNoise);
hold on;
xlabel('Range Bin');
ylabel('Amplitude [V]');
title('Typical Gaussian Noise Signal');