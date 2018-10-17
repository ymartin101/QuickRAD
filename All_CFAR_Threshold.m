%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% All_CFAR_Threshold.m: Simulate four CFAR detectors with Gaussian data
%% CA-, GOCA-, OS- and OSGO-CFAR detectors
clc;
close all;
clear;

% define inputs and parameters
D = 500;                                % length of data; number of samples
N = 16;                                 % length of window
Np = 1;                                 % number of pulses
j = 1i;                                 % use j as sqrt(-1)
PFA = 10^-6;                            % desired PFA
SNRdB = 20;                             % SNR in Decibels
SNR = 10^(SNRdB/10);                    % Decibels to linear SNR
iterations = 1000;
TCA_array = zeros(iterations,D);
TGOCA_array = zeros(iterations,D);
TOS_array = zeros(iterations,D);
TOSGO_array = zeros(iterations,D);
signal_array = zeros(iterations,D);

% include multiple targets
t = 1;                                  % number of targets and interfering targets

% include clutter edge
v = 5;                                  % variance of noise in clutter edge; v = 1 => no clutter edge
d = -1;                                % number of samples from centre to clutter edge start (distance)

for iteration = 1:iterations
    % 1xD matrix of complex Gaussian noise: (I + jQ)/sqrt(2); v scales second part => clutter edge
    noise = [((randn(1,(D/2) + d) + j.*randn(1,(D/2) + d))/sqrt(2)),(sqrt(v)*(randn(1,(D/2) - d) + j.*randn(1,(D/2) - d))/sqrt(2))];

    target_voltage = sqrt(SNR);
    target_signal = target_voltage*(randn(1,t) + j.*randn(1,t))/sqrt(2);
    signal1 = noise;	% noise only; H0 voltage
    for y = 0:(t - 1)   % update signal with targets from centre of data at intervals of 3 samples
        signal1((D/2) + 1 + 3*y) = signal1((D/2) + 1 + 3*y) + target_signal(y + 1);  % H1 voltage
    end

    % square law detector after target is added
    signal = abs(signal1).^2;

    %% CA-CFAR
    % compute CA-CFAR interference statistic 'g'
    g = zeros(D,1);
    for CUT = (N/2 + 1):1:(D - N/2)    % first N/2 and last N/2 samples in signal will NOT become the CUT (eg. 9 to 392)
        g(CUT) = mean([signal((CUT - N/2):(CUT - 1)) signal((CUT + 1):(CUT + N/2))]); % cell-average data in window (excluding CUT)
    end

    % compute CA-CFAR constant 'a' for desired PFA
    a = N.*(PFA.^(-(1./N)) - 1);

    % compute detection threshold 'T'
    T = a.*g;
    TCA_array(iteration,:) = T;

    %% GOCA-CFAR
    % compute GOCA-CFAR interference statistic 'g'
    g_GOCA = zeros(D,1);
    for CUT = (N/2 + 1):1:(D - N/2)    % first N/2 and last N/2 samples in signal will NOT become the CUT (eg. 9 to 392)
        Leading_window = sum(signal((CUT + 1):(CUT + N/2)));
        Lagging_window = sum(signal((CUT - N/2):(CUT - 1)));
        if Leading_window > Lagging_window
            g_GOCA(CUT) = Leading_window;
        else
            g_GOCA(CUT) = Lagging_window;
        end
    end

    % GOCA-CFAR constant 'a' for desired PFA
    a_GOCA = 2.4195;

    % compute detection threshold 'T'
    T_GOCA = a_GOCA.*g_GOCA;
    TGOCA_array(iteration,:) = T_GOCA;

    %% OS-CFAR
    % compute OS-CFAR interference statistic 'g'
    k_OS = 12;
    g_OS = zeros(D,1);
    for CUT = (N/2 + 1):1:(D - N/2)    % first N/2 and last N/2 samples in signal will NOT become the CUT (eg. 9 to 392)
        Window_sorted = sort([signal((CUT - N/2):(CUT - 1)) signal((CUT + 1):(CUT + N/2))]);    % sort N reference window samples in ascending order
        g_OS(CUT) = Window_sorted(k_OS);
    end

    % OS-CFAR constant 'a' for desired PFA
    a_OS = 20.9542;

    % compute detection threshold 'T'
    T_OS = a_OS.*g_OS;
    TOS_array(iteration,:) = T_OS;

    %% OSGO-CFAR
    % compute OSGO-CFAR interference statistic 'g'
    k_OSGO = 7;
    g_OSGO = zeros(D,1);
    for CUT = (N/2 + 1):1:(D - N/2)    % first N/2 and last N/2 samples in signal will NOT become the CUT (eg. 9 to 392)
        lagging_sorted = sort(signal((CUT - N/2):(CUT - 1)));    % sort lagging window samples in ascending order
        leading_sorted = sort(signal((CUT + 1):(CUT + N/2)));    % sort leading window samples in ascending order
        if leading_sorted(k_OSGO) > lagging_sorted(k_OSGO)
            g_OSGO(CUT) = leading_sorted(k_OSGO);
        else
            g_OSGO(CUT) = lagging_sorted(k_OSGO);
        end
    end

    % OSGO-CFAR constant 'a' for desired PFA
    a_OSGO = 12.85;

    % compute detection threshold 'T'
    T_OSGO = a_OSGO.*g_OSGO;
    TOSGO_array(iteration,:) = T_OSGO;
    
    % store signal in array
    signal_array(iteration,:) = signal;
end
    
% plot the threshold and noise
x = 0:(D - 1); % x-axis sample number
figure;
T = mean(10*log10(TCA_array),1);
T_GOCA =  mean(10*log10(TGOCA_array),1);
T_OS =  mean(10*log10(TOS_array),1);
T_OSGO =  mean(10*log10(TOSGO_array),1);
signal_mean =  mean(10*log10(signal_array),1);

plot(x,signal_mean,x,T,x,T_GOCA,x,T_OS,x,T_OSGO,'LineWidth',1);
hold on;
ylim([-60 70]);
xlabel('Range Bin');
ylabel('Amplitude [dB]');
title('Various CFAR Detectors');
lgd = legend('Data','CA Threshold','GOCA Threshold','OS Threshold','OSGO Threshold');
set(lgd,'Location','SouthEast');
