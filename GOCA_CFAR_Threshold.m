%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% GOCA_CFAR_Threshold.m: Simulate GOCA-CFAR detector with Gaussian data

%% GOCA-CFAR detector
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
v = 2;                                  % variance of noise in clutter edge; v = 1 => no clutter edge
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

% compute GOCA-CFAR interference statistic 'g'
g = zeros(D,1);
for CUT = (N/2 + 1):1:(D - N/2)    % first N/2 and last N/2 samples in signal will NOT become the CUT (eg. 9 to 392)
    Leading_window = sum(signal((CUT + 1):(CUT + N/2)));
    Lagging_window = sum(signal((CUT - N/2):(CUT - 1)));
    if Leading_window > Lagging_window
        g(CUT) = Leading_window;
    else
        g(CUT) = Lagging_window;
    end
end

% compute GOCA-CFAR constant 'a' for desired PFA
PFA_error = inf;

for a_i = 0:0.01:25
    PFA_summation = 0;
    for i = 0:(N/2 - 1)
        PFA_summation = PFA_summation + (factorial(N/2 - 1 + i)/(factorial(N/2 - 1).*factorial(i))).*(2 + a_i)^(-i);
    end
    PFA_i = 2*(((1 + a_i)^(-N/2)) -((2 + a_i)^(-N/2))*PFA_summation);
    if abs(PFA_i - PFA) < abs(PFA_error)
        a = a_i;
        PFA_error = abs(PFA_i - PFA);
    end
end

% compute detection threshold 'T'
T = a.*g;

% plot the threshold and noise
x = 0:(D - 1); % x-axis sample number
figure;
plot(x,10*log10(signal),x,10*log10(T));
xlabel('Range Bin');
ylabel('Amplitude [dB]');
ylim([-60 60]);
title('GOCA-CFAR Detector');
lgd = legend('Data','GOCA-CFAR Threshold Value');
set(lgd,'Location','SouthEast');
