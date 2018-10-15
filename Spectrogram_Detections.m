%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% Spectrogram_Detections.m: Generate spectrogram with detections

%% Startup
close all;
clc;
load('TFC15_008.mat');          % loads variables: Cdata, NumOfPRIs, NumOfRangeBins, PRI_s; must be included in .mat file

%% Parameters
c0 = 299792458;                             % speed of light
lambda = c0/(9*(10^9));                     % wavelength of the radar; f = 9 GHz
bandwidth = 10^7;                           % pulse width = 100us; pure sinusoid at 9 GHz
resolution_meters = c0/(2*bandwidth);       % range resolution in meters
j = sqrt(-1);                               % sqrt(-1)

% Sampling frequency, period
PRI = PRI_s;
PRF = 1/PRI;
ts = 1/PRF;             % unused? Equal to PRI?
target_doppler = PRF/4;	% all simulated target files use a target doppler of PRF/4

% Spectrogram
window_length = 512;    % length of window for FFT
overlap = 256;          % overlap in samples
NFFT = window_length;
freq_axis = (-(NFFT/2):1:((NFFT/2) - 1))*PRF/NFFT;

% OSGO-CFAR
a = 15;
N = 16;
k = round(5*N/12);

%% CFAR detection
RangeBin = 40;
X = Cdata(1:NumOfPRIs,RangeBin);
[S,F,T1,P] = spectrogram(X,window_length,overlap,NFFT,PRF); 
kc = floor((length(X) - overlap)/(window_length - overlap));
signal = fftshift(P);
signal2 = fftshift(P,1);

% matrix containing background statistics
g = zeros(NFFT,kc);                 % NFFT x kc matrix

% CUT and reference cells created + shifted for each time instance
for CUT = ((N/2) + 1):1:(NFFT - (N/2))        
    lagging_sorted = sort(signal(((CUT - N/2):(CUT - 1)),:));    % sort lagging window samples in ascending order
    leading_sorted = sort(signal(((CUT + 1):(CUT + N/2)),:));    % sort leading window samples in ascending order

    OSGO_diff = leading_sorted(k,:) - lagging_sorted(k,:);
    OSGO_lag_greater = find(OSGO_diff <= 0);
    OSGO_lead_greater = find(OSGO_diff > 0); 
    g(CUT,OSGO_lead_greater) = leading_sorted(k,OSGO_lead_greater);
    g(CUT,OSGO_lag_greater) = lagging_sorted(k,OSGO_lag_greater);
end

% Extract the data of the detected elements
T = a.*g;                                   % CFAR detection threshold; NFFT x kc matrix
signal_minus_T = signal - T;                % NFFT x kc matrix
detections_index = find(signal_minus_T > 0);% returns indices (counting down columns) of values > 0

detection_result = zeros(NFFT,kc);
detection_result(detections_index) = signal_minus_T(detections_index);  % same as signal_minus_T, but non-detections are all zero

%% Plot spectrogram and detections
figure;
imagesc(T1,freq_axis,10*log10(signal));
set(gca,'YDir','normal')    % correct orientation of image plot
colorbar;
colormap(jet);
hold on;
title('Clustered Detections in Measured Data');
xlabel('Time [seconds]');
ylabel('Frequency [Hz]');
x = overlap*PRI*(floor(detections_index/NFFT)+1);               % detections in time
y = floor(mod(detections_index,NFFT))*(PRF/NFFT) - (PRF/2);     % detections in frequency
plot(x,y,'k.');

