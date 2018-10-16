%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% OSGO_Find_a.m: Use OSGO-CFAR with desired PFA to find a

%% Startup
close all;
clc;
% Ensure "z" variable is loaded into workspace after generating the clutter data
PFA_desired = 10^-3;                        % desired PFA

%% Parameters
c0 = 299792458;                             % speed of light
lambda = c0/(6.9*(10^9));                   % wavelength of the radar; f = 6.9 GHz
bandwidth = 10^7;                           % pulse width = 100us; pure sinusoid at 6.9 GHz
resolution_meters = c0/(2*bandwidth);       % range resolution in meters
j = sqrt(-1);                               % sqrt(-1)

% Sampling frequency, period
PRF = 5000;             % assume PRF of 5000 Hz
ts = 1/PRF;             % unused? Equal to PRI?

% Spectrogram
window_length = 512;    % length of window for FFT
NFFT = window_length;
overlap = NFFT/2;       % overlap in samples
freq_axis = (-(NFFT/2):1:((NFFT/2) - 1))*PRF/NFFT;

% OSGO-CFAR parameters; tune until desired PFA is achieved
% "a" is unknown; determined iteratively later
N = 8;
k = round(5*N/12);
CT = 0.05425;

% Determine kc; frequency-domain signal will be of size NFFT x kc
kc = floor((length(z.') - overlap)/(window_length - overlap));
summation = zeros(kc,1);              % size of summation matrix

%% CFAR detection
[S,F,T1,P] = spectrogram(z.',window_length,overlap,NFFT,PRF); 
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
PFA_error = inf;        % initial error between actual and desired PFA
PFA_desired = 10^-3;    % desired PFA
summation = zeros(kc,1);% size of summation matrix

for a_i = 0:0.1:50     % test different "a" values to minimise PFA error    
    T = a_i.*g;                                 % CFAR detection threshold; NFFT x kc matrix
    signal_minus_T = signal - T;                % NFFT x kc matrix
    detections_index = find(signal_minus_T > 0);% returns indices of values > 0
    detection_result = zeros(NFFT,kc);
    detection_result(detections_index) = signal_minus_T(detections_index);
    summation = fftshift(sum(detection_result,1),1);	% summation for clustering
    detections_index = find(summation > CT);
    PFA_actual = length(detections_index)/length(summation);
    
    if abs(PFA_actual - PFA_desired) < PFA_error
        a = a_i;    % assign new alpha if PFA_actual has improved
        PFA_error = abs(PFA_actual - PFA_desired);
    end
end
