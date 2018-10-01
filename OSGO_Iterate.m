%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% OSGO_Iterate.m: Iterate using OSGO-CFAR to measured datasets for a varying parameter

%% FUNCTION: Used to determine PFA in measured target data for a varying parameter
function PFA = OSGO_Iterate(N)
%% Startup
close all;
clc;
load('TFC15_008.mat');	% input data with variables: Cdata, NumOfPRIs, NumOfRangeBins, PRI_s; must be included in .mat file

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
a = 20;
% N = 16;
k = round(5*N/12);
CT = 0.005;

% Determine kc; frequency-domain signal will be of size NFFT x kc
X = Cdata(1:NumOfPRIs,1);
kc = floor((length(X) - overlap)/(window_length - overlap));
summation = zeros(kc,NumOfRangeBins);   % size of summation matrix

%% CFAR detection
for RangeBin = 1:NumOfRangeBins
    X = Cdata(1:NumOfPRIs,RangeBin);
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

    % Extract the data of the detected elements
    T = a.*g;                                   % CFAR detection threshold; NFFT x kc matrix
    signal_minus_T = signal - T;                % NFFT x kc matrix
    detections_index = find(signal_minus_T > 0);% returns indices (counting down columns) of values > 0
    
    detection_result = zeros(NFFT,kc);
    detection_result(detections_index) = signal_minus_T(detections_index);  % same as signal_minus_T, but non-detections are all zero
    summation(:,(RangeBin)) = fftshift(sum(detection_result,1),1);          % sum the columns; produce row; becomes column in summation
    
    % Print progress percentage
    clc;
    fprintf('Processing: %0.2f %% complete \n',RangeBin*100/NumOfRangeBins);
end

%% Plot detections
range_bins = 1:1:size(Cdata,2);
time = (1:1:size(Cdata,1))*PRI;
% Plot data in time vs range bin
figure;
imagesc(range_bins,time,20*log10(abs(Cdata)));
colorbar;
colormap(jet);
hold on;
title('Detections in Measured Data');
xlabel('Range Bin');
ylabel('Time [seconds]');

% Apply clustering threshold and plot detection dots
detections_index = find(summation > CT);              % max number of entries = kc x NumOfRangeBins 
x = floor(detections_index./kc) + 1;                    % max x = NumOfRangeBins
y = PRI*overlap*floor(mod(detections_index,kc));        % max mod is (kc - 1) therefore add 1; max y = kc*PRI*overlap = max time
plot(x,y,'k.','MarkerSize',6);
colorbar;
colormap(jet);

% Calculate PFA for specific range bins
PFA = length(find(summation(:,[1:15,60:96]) > CT))/(length([1:15,60:96])*kc);

end