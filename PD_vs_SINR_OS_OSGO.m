%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% PD_vs_SINR_OS_OSGO.m: Generate PD vs SINR ROC for OS-, OSGO-CFAR

%% Setup and Inputs
clc;
close all;
clear;

% define inputs and parameters
D = 800;                                % length of data; number of samples
N = 16;                                 % length of window
j = 1i;                                 % use j as sqrt(-1)
PFA_desired = 10^-6;                    % desired PFA
k_factor = round(3*N/4);                % OS-CFAR k-factor
k_OSGO = round(5*N/12);
SNR_dB = 0:0.5:40;                      % SNR range in Decibels
SNR = 10.^(SNR_dB./10);                 % Decibels to linear SNR
k = 10^6;

%% OSGO-CFAR Setup and Constant

OSGO_a = 12.85;                        	% for PFA = 10^-6, N = 32
OSGO_PD_matrix = zeros(length(SNR),1);
OSGO_false_alarms = 0;

%% OS-CFAR Setup and Constant

OS_a = 20.9542;                         % for PFA = 10^-6, N = 32
OS_PD_matrix = OSGO_PD_matrix;
OS_false_alarms = 0;

%% Create noise

% complex Gaussian noise (I + jQ)/sqrt(2)
lagging_window = (randn(N/2,k) + j.*randn(N/2,k))/sqrt(2); % (N/2)xk

% complex Gaussian noise (I + jQ)/sqrt(2)
leading_window = (randn(N/2,k) + j.*randn(N/2,k))/sqrt(2); % (N/2)xk

%% Algorithm Thresholds

% OSGO threshold
% determine which ordered window's mean is larger, then add those cell entries into OSGO_g array
OSGO_lag_sorted = sort(real(lagging_window).^2 + imag(lagging_window).^2,1,'ascend');
OSGO_lead_sorted = sort(real(leading_window).^2 + imag(leading_window).^2,1,'ascend');
OSGO_diff = OSGO_lead_sorted(k_OSGO,:) - OSGO_lag_sorted(k_OSGO,:);
OSGO_lag_greater = find(OSGO_diff < 0);
OSGO_lead_greater = find(OSGO_diff >= 0);
OSGO_g(OSGO_lead_greater) = OSGO_lead_sorted(k_OSGO,OSGO_lead_greater);
OSGO_g(OSGO_lag_greater) = OSGO_lag_sorted(k_OSGO,OSGO_lag_greater);

OSGO_T = OSGO_a.*OSGO_g;

% OS threshold    
ref_cells = [lagging_window.', leading_window.'].';     % transpose and add windows into one matrix
OS_sorted = sort(real(ref_cells).^2 + imag(ref_cells).^2,1,'ascend');
OS_g = OS_sorted(k_factor,:);                           % assign g based on OS k-factor

OS_T = OS_a.*OS_g;

%% Algorithm PFAs
% H0 hypothesis (noise only)
H0 = (randn(1,k) + j*randn(1,k))/sqrt(2);

% OSGO false alarms
new_false_alarms = length(find(((real(H0).^2 + imag(H0).^2) - OSGO_T) > 0));
OSGO_false_alarms = OSGO_false_alarms + new_false_alarms;   % add up false alarms across all iterations

% OS false alarms
new_false_alarms = length(find(((real(H0).^2 + imag(H0).^2) - OS_T) > 0));
OS_false_alarms = OS_false_alarms + new_false_alarms;   % add up false alarms across all iterations

%% Algorithm PDs

for snr = 1:length(SNR)
    % H1 hypothesis (target and noise)
    H1 = H0 + ((randn(1,k) + j.*randn(1,k))/sqrt(2)).*sqrt(SNR(snr));

% OSGO detections
    new_detections = find(((real(H1).^2 + imag(H1).^2) - OSGO_T) > 0);
    OSGO_PD_matrix(snr,1) = length(new_detections)./k;   % add new PD to PD_matrix

% OS detections
    new_detections = find(((real(H1).^2 + imag(H1).^2) - OS_T) > 0);
    OS_PD_matrix(snr,1) = length(new_detections)./k;   % add new PD to PD_matrix
end

%% Plots and Results    
figure;
title('PD vs SINR');
xlabel('SINR [dB]');
ylabel('Probability of Detection');
grid on;
hold on;
plot(SNR_dB,OSGO_PD_matrix,'-.',SNR_dB,OS_PD_matrix,'--','LineWidth',1.5);
SNR_vs_PD_legend = legend('OSGO-CFAR','OS-CFAR');
set(SNR_vs_PD_legend,'Location','SouthEast');

% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 5 2.5];
% print('-depsc','-r0','C:\Users\mrtmo\Desktop\UCT\EEE4022S\Report\Images\OSGO_roc')
