%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% Target_Gaussian.m: Generate simulated target in Gaussian noise

%% Produce a target in Gaussian noise
% Produce Gaussian noise
D = 500;                                % Data set size
I = randn(1,D);                         % In-phase component
Q = randn(1,D);                         % Quadrature component
GaussianNoise = (I + 1i.*Q)/sqrt(2);	% Produce 1xD matrix of Gaussian noise

% Produce a Gaussian target
SNR_dB = 20;							% SNR in dB
SNR = 10^(SNR_dB/10);					% SNR linear
Target = sqrt(SNR).*(randn(1,1) + 1i.*randn(1,1))/sqrt(2);	% Target data

% Insert target at centre of data
signal = GaussianNoise;
signal(D/2) = signal(D/2) + Target;
signal = signal.^2;

% Plot data (dB)
figure;
n = 0:(D - 1);
plot(n,20*log10(signal));
hold on;
xlabel('Range Bin');
ylabel('Amplitude [dB]');
title('Target in Noise Signal');
SNR_dB = 20;							% SNR in dB
SNR = 10^(SNR_dB/10);                   % SNR linear
Target = sqrt(SNR).*(randn(1,1) + 1i.*randn(1,1))/sqrt(2);	% Target data