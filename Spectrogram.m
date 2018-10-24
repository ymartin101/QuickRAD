%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% Spectrogram.m: Generate spectrogram and inserted a simulated target in one range bin

%% Setup
clear all;
close all;
% select input dataset
load('TFC15_008');  % loads a variables Cdata, Info, NumOfPRIs, NumOfRangeBins, PCI, PRI_s, Waveform, lamda  
RangeBin = 90;      % select range bin to process

%% Parameters and calculations
PRI = PRI_s; % from dataset
PRF = 1/PRI;
ts = 1/PRF;
num_range_lines = size(Cdata,1);
start_range = num_range_lines/2;
start_time = start_range*PRI;
stop_range = num_range_lines/4*3;
stop_time = stop_range*PRI;
H0 = Cdata(start_range:stop_range,RangeBin);
H1 = H0;
num_samples = length(H0);

%% Insert target (comment out if you want the spectrogram for a measured target in a specific range bin)

  SINR = 0; % SINR in dB
  target_doppler = 600;  % specify between -PRF/2 to PRF/2 

  % Estimate power of noise and clutter
  H0_power = var(H0); % assume Gaussian PDF
  target_power = (H0_power)*10^(SINR/10); % target power to linear 
  
  % Steering vector to incorporate speed of the target
  m = (0:1:(num_samples-1)).';
  vd = target_doppler/PRF; 
  v = exp(1j*2*pi*m*vd); 
  
  % Target signal
  target_amplitude = sqrt(target_power/2)*(randn(1,1)+1j*randn(1,1)); % target amplitude
  target_signal = repmat(target_amplitude,num_samples,1).*v;            % target signal
  
  % Add the simulated target to the measured clutter
  H1 = H0 + target_signal;
  
%% Generate and plot Spectrogram
window_length = 512;
overlap = 256;          % number of overlapping window samples
NFFT = window_length;
freq_axis = (-NFFT/2:1:(NFFT/2-1))*PRF/NFFT;
[S,F,T1,P] = spectrogram(H1,window_length,overlap,NFFT,PRF);  % S = SPECTROGRAM(H0Data,WINDOW,NOVERLAP,NFFT) 

figure; 
P_fftshift = fftshift(P, 1);
imagesc((T1+start_time),freq_axis,10*log10(P_fftshift)); 
axis xy
xlabel('Time [seconds]'); ylabel('Frequency [Hz]');
title('Spectrogram Plot');
colorbar;
colormap(jet);
