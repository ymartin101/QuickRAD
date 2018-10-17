%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% PD_Extended_OS.m: Determine PD for a varying parameter in a "extended" dataset using OS-CFAR

%% FUNCTION: Used to determine PFA and PD for a varying parameter
function PD_mean = PD_Extended_OS(target_doppler)
%% Add extended target to a clutter dataset and determine PD
clc;
close all;

%% Inputs
load('CFA17_003.mat');  % loads variables: Cdata, Info, NumOfPRIs, NumOfRangeBins, PCI, PRI_s, Waveform
a = 25;
N = 8;
k = round(3*N/4);
NFFT = 512;             % length of window for FFT
clustering_threshold = 0.05;
% target_doppler = 600;

%% Setup
% Sampling frequency, period
PRI = PRI_s;
PRF = 1/PRI;
target_length = 5;      % number of range bins for target to occupy
fD = ceil(((target_doppler) + (PRF/2))*NFFT/PRF);   % positive frequency (target coming towards radar)
fD1 = (NFFT/2) - (fD - (NFFT/2));   % negative frequency (target moving away from radar)
kc1_mat5 = zeros(NumOfRangeBins,1);
kc2_mat1 = zeros(NumOfRangeBins,1);
PD_mat1 = zeros(NumOfRangeBins,1);
PD_mat2 = zeros(NumOfRangeBins,1);
PD_mat3 = zeros(NumOfRangeBins,1);

%% Parameters
c0 = 299792458;             % speed of light
lambda = c0/(6.9*(10^9));	% wavelength of the radar; f = 6.9 GHz for the clutter datasets
v_target = lambda*target_doppler/2; % target speed [m/s]
dt = 15/v_target;           % v_target = dRangeBin/dt = 15/dt (1 RangeBin = 15 m)

% Spectrogram
overlap = NFFT/2;          % overlap in samples
window_length = NFFT;

% Get number of target cells
X = Cdata(1:NumOfPRIs,1);
kc = floor((length(X) - overlap)/(window_length - overlap));

% For target insertion
C_backup = Cdata;
t1 = 0;
t2 = 0;
SINR = 0;   % dB
break_loop = 0;
RBMax = 0;

%% Range bins with only clutter
if abs(target_doppler) <= 20
    RangeBin = 1:5;
    RBMax = 5;
    t1 = 1;
    t2 = NumOfPRIs;	% update t2 with new t1 (PRIs)
    kc1 = t1*kc/NumOfPRIs;
    kc1_mat5(RangeBin) = kc1;
    kc2 = t2*kc/NumOfPRIs;
    kc2_mat1(RangeBin) = kc2;
    
    % Target 
    H0 = Cdata(:,RangeBin);	% column vector from Cdata; use RangeBin
    num_samples = length(H0);
    
    % Estimate power of noise and clutter
    H0_power = var(H0);           % assume Gaussian PDF (incorrect, but just a first order estimate) 

    % Compute target power
    target_power = (H0_power)*10^(SINR/10); % Target power - linear 

    % Steering vector is used to incorporate the speed of the target
    m = (0:1:(num_samples - 1)).';
    vd = -target_doppler/PRF; 
    v = exp(1j*2*pi*m*vd); 

    % Calculate target signal
    target_amplitude = sqrt(target_power);       	% target amplitude (Swerling I model)
    target_signal = repmat(target_amplitude,num_samples,1).*v;                      % target signal (Swerling I model)

    % Add the simulated target to the measured clutter
    H1 = H0 + target_signal;
    C_backup(:,RangeBin) = H1;
    C_edited(:,RangeBin) = 1;
    
else
    % insert target into five range bins
    for rbin = 1:target_length  % 1:length(desired length of target)
        t1 = 0;
        t2 = 0;
        break_loop = 0;
        for RangeBin = rbin:NumOfRangeBins
            t1 = ceil((t2/PRF)*PRF);	% update t1 with previous t2 (PRIs)
            t2 = (t1/PRF + dt)*PRF;     % gradient of v_target
            if t2 > size(Cdata,1)
                t2 = NumOfPRIs;
                break_loop = 1;
            end
            % transform time values to kc values
            if rbin == 5    % last range bin portion of extended target
                kc1 = t1*kc/NumOfPRIs;
                kc1_mat5(RangeBin) = kc1;
                RBMax = RangeBin;
            end
            if rbin == 1    % first range bin portion of extended target
                kc2 = t2*kc/NumOfPRIs;
                kc2_mat1(RangeBin) = kc2;
            end

            H0 = Cdata(t1:t2,RangeBin);	% column vector from Cdata; use RangeBin
            num_samples = length(H0);
            H0_power = var(H0);           % assume Gaussian PDF (incorrect, but just a first order estimate) 

            % Compute target power
            target_power = (H0_power)*10^(SINR/10); % target power - linear 

            % Steering vector is used to incorporate the speed of the target
            m = (0:1:(num_samples - 1)).';
            vd = -target_doppler/PRF; 
            v = exp(1j*2*pi*m*vd); 

            % Calculate target signal
            target_amplitude = sqrt(target_power);	% target amplitude (Swerling I model)
            target_signal = repmat(target_amplitude,num_samples,1).*v;	% target signal (Swerling I model)

            % Add the simulated target to the measured clutter
            H1 = H0 + target_signal;
            C_backup(t1:t2,RangeBin) = H1;
            if break_loop == 1
                break
            end
        end
    end
end

%% OSGO-CFAR on clutter only; for PD
% CFAR detection
for RangeBin = 1:NumOfRangeBins
    X = C_backup(1:NumOfPRIs,RangeBin);
    [S,F,T1,P] = spectrogram(X,window_length,overlap,NFFT,PRF); 
    signal = fftshift(P,1);
    
    % matrix containing background statistics
    g = zeros(NFFT,kc);	% NFFT x kc matrix

    % CUT and reference cells created + shifted for each time instance
    for CUT = ((N/2) + 1):1:(NFFT - (N/2))        
        window_sorted = sort([signal(((CUT + 1):(CUT + N/2)),:); signal(((CUT - N/2):(CUT - 1)),:)]); % sort window samples in ascending order
        g(CUT,:) = window_sorted(k,:);
    end

    % Extract the data of the detected elements
    T = a.*g;                                   % CFAR detection threshold; NFFT x kc matrix
    signal_minus_T = signal - T;                % NFFT x kc matrix
    detections_index = find(signal_minus_T > 0);% returns indices (counting down columns) of values > 0
    
    %% Determine PD; consider all cases of target reaching last range bin, etc. (depends on speed)
    if RangeBin <= RBMax
        if RangeBin > target_length && fD ~= 128    % find PD; length of target = 5
            difference1 = signal_minus_T(fD1+1,kc1_mat5(RangeBin):kc2_mat1(RangeBin));
            if fD < NFFT
                difference2 = signal_minus_T(fD1,kc1_mat5(RangeBin):kc2_mat1(RangeBin));
                difference3 = signal_minus_T(fD1-1,kc1_mat5(RangeBin):kc2_mat1(RangeBin));
                diff1_max = max([length(find(difference1 > 0)),length(find(difference2 > 0)),length(find(difference3 > 0))]);
            else
                diff1_max = length(find(difference1 > 0));
            end
            target_detections = diff1_max;
            PD_mat1(RangeBin) = target_detections/length(kc1_mat5(RangeBin):kc2_mat1(RangeBin));
        elseif RangeBin <= target_length
            difference1 = signal_minus_T(fD1+1,1:kc2_mat1(RangeBin));
            if fD < NFFT
                difference2 = signal_minus_T(fD1,1:kc2_mat1(RangeBin));
                difference3 = signal_minus_T(fD1+2,1:kc2_mat1(RangeBin));
                diff2_max = max([length(find(difference1 > 0)),length(find(difference2 > 0)),length(find(difference3 > 0))]);
            else
                diff2_max = length(find(difference1 > 0));
            end
            target_detections = diff2_max;
            PD_mat3(RangeBin) = target_detections/length(1:kc2_mat1(RangeBin));
        end
    end
    
    clc;
    fprintf('Processing: %0.2f %% complete \n',RangeBin*100/NumOfRangeBins);
end

PD_mean1 = mean(PD_mat1(PD_mat1 > 0));
PD_mean2 = mean(PD_mat2(PD_mat2 > 0));
PD_mean3 = mean(PD_mat3(PD_mat3 > 0));
PD_mean = max([PD_mean1,PD_mean2,PD_mean3]);    % average PD from all range bins

end