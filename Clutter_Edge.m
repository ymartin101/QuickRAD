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

% 1xD matrix of complex Gaussian noise: (I + jQ)/sqrt(2); v scales second part => clutter edge
noise = [((randn(1,(D/2) + d) + j.*randn(1,(D/2) + d))/sqrt(2)),(sqrt(v)*(randn(1,(D/2) - d) + j.*randn(1,(D/2) - d))/sqrt(2))];

% insert target(s)
signal1 = noise; % noise only; H0

% square law detector after target is added
signal = abs(signal1).^2;

% plot the relevant information
x = 0:(D - 1); % x-axis sample number
figure;
plot(x,10*log10(signal));
hold on;
xlabel('Range Bin');
ylabel('Amplitude [dB]');
title('Typical Clutter Edge');
