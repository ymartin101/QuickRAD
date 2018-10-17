%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% OSGO_CFAR_Threshold.m: Simulate OSGO-CFAR detector with Gaussian data

%% OSGO-CFAR detector
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
k = round(5*N/12);                      % k-factor for OS-CFAR

% include multiple targets
t = 2;                                  % number of targets and interfering targets (max 3); default is 1

% include clutter edge
v = 1;                                  % variance of noise in clutter edge; v = 1 => no clutter edge (default)
d = -1;                                  % number of samples from centre to clutter edge start (distance)

% 1xD matrix of complex Gaussian noise: (I + jQ)/sqrt(2); v scales second part => clutter edge
noise = [((randn(1,(D/2) + d) + j.*randn(1,(D/2) + d))/sqrt(2)),(sqrt(v)*(randn(1,(D/2) - d) + j.*randn(1,(D/2) - d))/sqrt(2))];
    
target_voltage = sqrt(SNR);
target_signal = target_voltage*(randn(1,t) + j.*randn(1,t))/sqrt(2);
signal1 = noise; % noise only; H0
for y = 0:(t - 1)   % update signal with targets from centre of data at intervals of 3 samples
    signal1((D/2) + 1 + 3*y) = signal1((D/2) + 1 + 3*y) + target_signal(y + 1);  % H1
end

% square law detector after target is added
signal = abs(signal1).^2;

% compute OSGO-CFAR interference statistic 'g'
g = zeros(D,1);
for CUT = (N/2 + 1):1:(D - N/2)    % first N/2 and last N/2 samples in signal will NOT become the CUT (eg. 9 to 392)
    lagging_sorted = sort(signal((CUT - N/2):(CUT - 1)));    % sort lagging window samples in ascending order
    leading_sorted = sort(signal((CUT + 1):(CUT + N/2)));    % sort leading window samples in ascending order
    if leading_sorted(k) > lagging_sorted(k)
        g(CUT) = leading_sorted(k);
    else
        g(CUT) = lagging_sorted(k);
    end
end

% compute OS-CFAR constant 'a' for desired PFA
PFA_error = inf;

for a_i = 0:0.05:100
    PFA_temp = 0;
    
    for m = 0:((N/2) - k)
        for p = 0:((N/2) - k)
            PFA_temp = PFA_temp + (factorial(((N/2) - k))./(factorial(m).*factorial((N/2) - k - m))).*...
                (factorial(((N/2) - k))./(factorial(p).*factorial((N/2) - k - p))).*...
                ((-1).^(N - (2*k) - m - p)./((N/2) - p)).*...
                (gamma(N - m - p).*gamma(a_i + 1)./gamma(N - m - p + a_i + 1));
        end
    end
    
    PFA_temp = (2*(k.^2)*(factorial(N/2)./(factorial((N/2) - k).*factorial(k))).^2).*PFA_temp;
    
    if abs(PFA_temp - PFA) < PFA_error
        a = a_i;
        PFA_error = abs(PFA_temp - PFA);
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
ylim([-80 80]);
title('OSGO-CFAR Detector');
lgd = legend('Data','OSGO-CFAR Threshold Value');
set(lgd,'Location','SouthEast');
