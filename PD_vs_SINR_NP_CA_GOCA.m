%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% PD_vs_SINR_NP_CA_GOCA.m: Generate PD vs SINR ROC for CA-, GOCA-CFAR and NP

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
SNR_dB = 0:0.5:30;                      % SNR range in Decibels
SNR = 10.^(SNR_dB./10);                   % Decibels to linear SNR
iterations = 10^6;                      % number of iterations/CUTs
k = 10^6;

%% Neyman-Pearson Setup
NP_false_alarms = 0;
NP_PD_matrix = zeros(length(SNR), ceil(iterations/10^6));

%% CA-CFAR Setup and Constant
CA_a = N*(PFA_desired^(-(1/N))-1);      % CA detection threshold
CA_false_alarms = 0;
CA_PD_matrix = NP_PD_matrix;

%% GOCA-CFAR Setup and Constant
GOCA_a = 2.4195;                        % for PFA = 10^-6, N = 16
GOCA_PD_matrix = NP_PD_matrix;
GOCA_false_alarms = 0;

%% Create noise
% complex Gaussian noise (I + jQ)/sqrt(2)
lagging_window = (randn(N/2,k) + j.*randn(N/2,k))/sqrt(2); % (N/2)xk

% complex Gaussian noise (I + jQ)/sqrt(2)
leading_window = (randn(N/2,k) + j.*randn(N/2,k))/sqrt(2); % (N/2)xk

%% Mean noise in reference cells
lagging_mean = mean((real(lagging_window).^2 + imag(lagging_window).^2),1);
leading_mean = mean((real(leading_window).^2 + imag(leading_window).^2),1);

%% Algorithm Thresholds 
% NP threshold
NP_T = -log(PFA_desired);

% CA-CFAR threshold
CA_g = (lagging_mean + leading_mean)/2;
CA_T = CA_a.*CA_g;

% GOCA threshold
sum_lagging = sum((real(lagging_window).^2 + imag(lagging_window).^2),1);
sum_leading = sum((real(leading_window).^2 + imag(leading_window).^2),1);
GOCA_g = zeros(1,k);

% determine which window's sum is greater, then add those cell entries into GOCA_g array
window_diff = sum_leading - sum_lagging;
lag_greater = find(window_diff <= 0);
lead_greater = find(window_diff > 0); 
GOCA_g(lead_greater) = sum_leading(lead_greater);
GOCA_g(lag_greater) = sum_lagging(lag_greater);

GOCA_T = GOCA_a.*GOCA_g;

%% H0 hypothesis (noise only)
H0 = (randn(1,k) + j*randn(1,k))/sqrt(2);

%% Algorithm PDs
for snr = 1:length(SNR)
    % H1 hypothesis (target and noise)
    H1 = H0 + ((randn(1,k) + j.*randn(1,k))/sqrt(2)).*sqrt(SNR(snr));

    % NP detections
    new_detections = find(((real(H1).^2 + imag(H1).^2) - NP_T) > 0);
    NP_PD_matrix(snr,1) = length(new_detections)./k;   % add new PD to PD_matrix

    % CA detections
    new_detections = find(((real(H1).^2 + imag(H1).^2) - CA_T) > 0);
    CA_PD_matrix(snr,1) = length(new_detections)./k;   % add new PD to PD_matrix

    % GOCA detections
    new_detections = find(((real(H1).^2 + imag(H1).^2) - GOCA_T) > 0);
    GOCA_PD_matrix(snr,1) = length(new_detections)./k;   % add new PD to PD_matrix
end

%% Plots and Results
figure;
title('PD vs SINR');
xlabel('SINR [dB]');
ylabel('Probability of Detection');
grid on;
hold on;
plot(SNR_dB,NP_PD_matrix,'-',SNR_dB,CA_PD_matrix,'--r',SNR_dB,GOCA_PD_matrix,'-.k','LineWidth',1.5);
SNR_vs_PD_legend = legend('NP','CA-CFAR','GOCA-CFAR');
set(SNR_vs_PD_legend,'Location','SouthEast');
