%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% Clutter_Edge.m: Generate simulated clutter edge in Gaussian noise

%% Setup
clc;
close all;
clear;

%% Define inputs and parameters
D = 800;                                % length of data; number of samples
N = 24;                                 % length of window
j = 1i;                                 % use j as sqrt(-1)
PFA = 10^-3;                            % desired PFA
SNR_dB = 20;                            % desired SNR in dB
SNR = 10^(SNR_dB/10);                   % convert dB to linear
t = 0;                                  % targets + interfering targets (max 3)
v = 10;                                 % variance of noise in clutter edge; v = 1 => no clutter edge
d = 0;                                  % number of samples from centre to clutter edge start (distance)

% 1xD matrix of complex Gaussian noise: (I + jQ)/sqrt(2); v scales second part => clutter edge
noise = [((randn(1,(D/2) + d) + j.*randn(1,(D/2) + d))/sqrt(2)),(sqrt(v)*(randn(1,(D/2) - d) + j.*randn(1,(D/2) - d))/sqrt(2))];

% mean noise power for square law
P_Noise = (mean(real(noise).^2) + mean(imag(noise).^2));

% noise variance before square law detector
V_noise = var(real(noise) + imag(noise));

% insert target(s)
target_voltage = sqrt(SNR);
target_signal = target_voltage*randn(1,t)*(1 + j)/sqrt(2);
signal = noise; % noise only; H0
for y = 0:(t - 1)   % update signal with targets from centre of data at intervals of 3 samples ???
    signal((D/2) + 1 + 3*y) = (signal((D/2) + 1 + 3*y) + target_signal(y + 1)).^2;  % H1
end

% compute CA-CFAR interference statistic 'g'
g = zeros(D,1);
for CUT = (N/2 + 1):1:(D - N/2)    % first N/2 and last N/2 samples in signal will NOT become the CUT (eg. 9 to 392)
    g(CUT) = mean([signal((CUT - N/2):(CUT - 1)) signal((CUT + 1):(CUT + N/2))]); % cell-average data in window (excluding CUT)
end

% compute CA-CFAR constant 'a' for desired PFA
a = N.*(PFA.^(-(1./N)) - 1);    % 1xk size

% compute detection threshold 'T'
T = a.*g;  %1xk

% plot the relevant information
x = 0:(D - 1); % x-axis sample number
figure;
plot(x,20*log10(signal));
hold on;
xlabel('Range Bin');
ylabel('Amplitude [dB]');
title('Typical Clutter Edge');
