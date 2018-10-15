%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% NP_Threshold.m: Neyman-Pearson detector for simulated data

%% Setup
clc;
close all;
clear;

% Define inputs and parameters
N = 200;                                % number of samples
j = 1i;                                 % use j as sqrt(-1)
PFA_desired = 10^-6;                    % desired PFA
SNR_dB = 20;                            % SNR in Decibels
SNR = 10^(SNR_dB/10);                   % Decibels to linear SNR
v = 1;                                  % default variance is 1; v = B^2
n = 0:(N - 1);

% complex Gaussian noise (I + jQ)/sqrt(2)
noise = (randn(1,N) + j.*randn(1,N))/sqrt(2);

% NP threshold
T = -log(PFA_desired)*sqrt(v);          % v = 1 by default

% constant threshold value line
x = T*ones(1,length(n));

% define target
target = sqrt(SNR)*(randn(1,1) + j.*randn(1,1))/sqrt(2);
signal = noise;

% insert target at centre of data
signal(N/2) = noise(N/2) + target;

% square law detector after target is added
signal = abs(signal).^2;

% plot data and NP threshold
figure;
plot(n,20*log10(signal),n,20*log10(x));
ylim([-60 60]);
hold on;
xlabel('Range Bin');
ylabel('Amplitude [dB]');
title('Neyman-Pearson Detector');
lgd = legend('Data','NP Threshold Value');
set(lgd,'location','NorthWest');
