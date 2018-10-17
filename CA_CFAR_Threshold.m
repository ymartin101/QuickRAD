%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% CA_CFAR_Threshold.m: Simulate CA-CFAR detector with Gaussian data

%% CA-CFAR detector
clc;
close all;
clear;

% define inputs and parameters
D = 250;                                % length of data; number of samples
N = 16;                                 % length of window
Np = 1;                                 % number of pulses
j = 1i;                                 % use j as sqrt(-1)
PFA = 10^-6;                            % desired PFA
SNRdB = 20;                             % SNR in Decibels
SNR = 10^(SNRdB/10);                    % Decibels to linear SNR

% include multiple targets
t = 1;                                  % number of targets and interfering targets (max 3)

% include clutter edge
v = 1;                                  % variance of noise in clutter edge; v = 1 => no clutter edge
d = -1;                                  % number of samples from centre to clutter edge start (distance)

% 1xD matrix of complex Gaussian noise: (I + jQ)/sqrt(2); v scales second part => clutter edge
noise = [((randn(1,(D/2) + d) + j.*randn(1,(D/2) + d))/sqrt(2)),(sqrt(v)*(randn(1,(D/2) - d) + j.*randn(1,(D/2) - d))/sqrt(2))];
    
% insert target(s)
target_voltage = sqrt(SNR);
target_signal = target_voltage*(randn(1,t) + j.*randn(1,t))/sqrt(2);
signal1 = noise; % noise only; H0
for y = 0:(t - 1)   % update signal with targets from centre of data at intervals of 3 samples
    signal1((D/2) + 1 + 3*y) = signal1((D/2) + 1 + 3*y) + target_signal(y + 1);  % H1
end

% square law detector after target is added
signal = abs(signal1).^2;

% compute CA-CFAR interference statistic 'g'
g = zeros(D,1);
for CUT = (N/2 + 1):1:(D - N/2)    % first N/2 and last N/2 samples in signal will NOT become the CUT (eg. 9 to 392)
    g(CUT) = mean([signal((CUT - N/2):(CUT - 1)) signal((CUT + 1):(CUT + N/2))]); % cell-average data in window (excluding CUT)
end

% compute CA-CFAR constant 'a' for desired PFA
a = N.*(PFA.^(-(1./N)) - 1);

% compute detection threshold 'T'
T = a.*g;

% plot the threshold and noise
x = 0:(D - 1); % x-axis sample number
figure;
plot(x,10*log10(signal),x,10*log10(T));
xlabel('Range Bin');
ylabel('Amplitude [dB]');
ylim([-60 60]);
title('CA-CFAR Detector');
lgd = legend('Data','CA-CFAR Threshold Value');
set(lgd,'Location','SouthEast');

