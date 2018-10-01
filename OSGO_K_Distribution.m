%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% OSGO_Find_a.m: Use OSGO-CFAR with desired PFA to find a

%% Startup
close all;
clc;
load('CFA17_002.mat');	% input data with variables: Cdata, NumOfPRIs, NumOfRangeBins, PRI_s; must be included in .mat file
RangeBin = 1;           % select range bin to process
PFA_desired = 10^-6;    % user-selected desired PFA

%% Parameters
c0 = 299792458;                             % speed of light
lambda = c0/(6.9*(10^9));                   % wavelength of the radar; f = 6.9 GHz
bandwidth = 10^7;                           % pulse width = 100us; pure sinusoid at 6.9 GHz
resolution_meters = c0/(2*bandwidth);       % range resolution in meters
j = sqrt(-1);                               % sqrt(-1)

% Sampling frequency, period
PRI = PRI_s;
PRF = 1/PRI;
ts = 1/PRF;             % unused? Equal to PRI?

% Spectrogram
window_length = 512;    % length of window for FFT
NFFT = window_length;
overlap = NFFT/2;       % overlap in samples
freq_axis = (-(NFFT/2):1:((NFFT/2) - 1))*PRF/NFFT;

% OSGO-CFAR parameters
% "a" is unknown
N = 8;
k = round(5*N/12);
CT = 0.01;

% Determine kc; frequency-domain signal will be of size NFFT x kc
X = Cdata(:);
kc = floor((length(X) - overlap)/(window_length - overlap));
summation = zeros(kc,1);              % size of summation matrix

%% CFAR detection
X = Cdata(:);
[S,F,T1,P] = spectrogram(X,window_length,overlap,NFFT,PRF); 
kc = floor((length(X) - overlap)/(window_length - overlap));
signal = fftshift(P,1);             % centre zero-frequency component  

% Matrix of background statistic
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

%% Calculate PFA and minimise error between actual and desired PFA
PFA_error = inf;
PFA_temp = 0;

for a_i = 40:0.5:80    % test different "a" values to minimise PFA error    
    T = a_i.*g;                                 % CFAR detection threshold; NFFT x kc matrix
    signal_minus_T = signal - T;                % NFFT x kc matrix
    detections_index = find(signal_minus_T > 0);% returns indices (counting down columns) of values > 0
    detection_result = zeros(NFFT,kc);
    detection_result(detections_index) = signal_minus_T(detections_index);
    summation = fftshift(sum(detection_result,1),1);	% summation for clustering
    detections_index = find(summation > CT);	% max number of entries = kc x NumOfRangeBins 
    PFA_temp = length(detections_index)/length(summation);
    
    if abs(PFA_temp - PFA_desired) < PFA_error
        a = a_i;
        PFA_error = abs(PFA_temp - PFA_desired);
    end
end

