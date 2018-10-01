%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% PD_Multiple_OS.m: Determine PD for a varying parameter in a "multiple" dataset using OS-CFAR

%% FUNCTION: Used to determine PFA and PD for a varying parameter
function PD_mean = PD_Multiple_OS(target_doppler)
%% Add multiple targets to a clutter dataset and determine PD
clc;
close all;

%% Inputs
load('CFA17_001.mat');  % loads variables: Cdata, Info, NumOfPRIs, NumOfRangeBins, PCI, PRI_s, Waveform
a = 25;
N = 16;
k = round(3*N/4);
NFFT = 512;             % length of window for FFT
clustering_threshold = 0.001;
% target_doppler = 600;

%% Setup
% Sampling frequency, period
PRI = PRI_s;
PRF = 1/PRI;
fD = ceil(((target_doppler) + (PRF/2))*NFFT/PRF);
fD1 = (NFFT/2) - (fD - (NFFT/2));
kc1_mat = zeros(NumOfRangeBins,1);
kc2_mat = zeros(NumOfRangeBins,1);
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
RBMax = -1;

% %% CFAR for PFA
% summation = zeros(kc,NumOfRangeBins);
% for RangeBin = 1:NumOfRangeBins
%     X = Cdata(1:NumOfPRIs,RangeBin);
%     [S,F,T1,P] = spectrogram(X,window_length,overlap,NFFT,PRF); 
%     kc = floor((length(X) - overlap)/(window_length - overlap));
%     signal = fftshift(P,1);
%     
%     % matrix containing background statistics
%     g = zeros(NFFT,kc);                 % NFFT x kc matrix
% 
%     % CUT and reference cells created + shifted for each time instance  
%     for CUT = ((N/2) + 1):1:(NFFT - (N/2))        
%         window_sorted = sort([signal(((CUT + 1):(CUT + N/2)),:); signal(((CUT - N/2):(CUT - 1)),:)]); % sort window samples in ascending order
%         g(CUT,:) = window_sorted(k,:);
%     end
% 
%     % Extract the data of the detected elements
%     T = a.*g;                                   % CFAR detection threshold; NFFT x kc matrix
%     signal_minus_T = signal - T;                % NFFT x kc matrix
%     detections_index = find(signal_minus_T > 0);% returns indices (counting down columns) of values > 0
%     
%     detection_result = zeros(NFFT,kc);
%     detection_result(detections_index) = signal_minus_T(detections_index);  % same as signal_minus_T, but non-detections are all zero
%     summation(:,(RangeBin)) = fftshift(sum(detection_result,1),1);          % sum the columns; produce row; becomes column in summation
%     
%     clc;
%     fprintf('Processing: %0.2f %% complete \n',RangeBin*100/NumOfRangeBins);
% end
% 
% % Calculate PFA
% detections_index = find(summation > clustering_threshold);          % max number of entries = kc x NumOfRangeBins 
% PFA = length(detections_index)/(kc*NumOfRangeBins);	% calculate probability of false alarm
% 
%% CFAR for PD
%% Range bins with only clutter
if abs(target_doppler) <= 20
    RangeBin = 1;
    t1 = 1;
    t2 = NumOfPRIs;	% update t2 with new t1 (PRIs)
    kc1 = t1*kc/NumOfPRIs;
    kc1_mat(RangeBin) = (kc1);
    kc1_mat(49 - RangeBin) = (kc1);
    kc2 = t2*kc/NumOfPRIs;
    kc2_mat(RangeBin) = (kc2);
    kc2_mat(49 - RangeBin) = (kc2);

%% Insert first target (towards radar)
    H0 = Cdata(:,(49 - RangeBin));	% column vector from Cdata; flipped plot
    num_samples = length(H0);
    H0_power = var(H0);             % assume Gaussian PDF (incorrect, but just a first order estimate) 

    % Compute target power
    target_power = (H0_power)*10^(SINR/10); % target power - linear 

    % Steering vector is used to incorporate the speed of the target
    m = (0:1:(num_samples - 1)).';
    vd = target_doppler/PRF; 
    v = exp(1j*2*pi*m*vd); 

    % Calculate target signal
    target_amplitude = sqrt(target_power);	% target amplitude (Swerling I model)
    target_signal = repmat(target_amplitude,num_samples,1).*v;	% target signal (Swerling I model)

    % Add the simulated target to the measured clutter
    H1 = H0 + target_signal;
    C_backup(:,(49 - RangeBin)) = H1;
    
%% Insert second target (away from radar)
    H0 = Cdata(:,RangeBin);	% column vector from Cdata; use RangeBin
    num_samples = length(H0);
    
    % Estimate power of noise and clutter
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
    C_backup(:,RangeBin) = H1;
    
else
    for RangeBin = 1:NumOfRangeBins
        t1 = ceil((t2/PRF)*PRF);   % update t1 with previous t2 (PRIs)
        t2 = (t1/PRF + dt)*PRF;   % gradient of v_target
        if t2 > size(Cdata,1)
            t2 = NumOfPRIs;
            RBMax = RangeBin;
            break_loop = 1;
        end
        kc1 = t1*kc/NumOfPRIs;
        if t1 == 0
            kc1 = 1;
        end
        % transform time values to kc values
        kc1_mat(RangeBin) = (kc1);
        kc2 = t2*kc/NumOfPRIs;
        kc2_mat(RangeBin) = (kc2);

    %% Insert first target (towards radar)
        H0 = Cdata(t1:t2,(49 - RangeBin));	% column vector from Cdata; flipped plot
        num_samples = length(H0);
        H0_power = var(H0); % assume Gaussian PDF (incorrect, but just a first order estimate) 

        % Compute target power
        target_power = (H0_power)*10^(SINR/10); % target power - linear 

        % Steering vector is used to incorporate the speed of the target
        m = (0:1:(num_samples - 1)).';
        vd = target_doppler/PRF; 
        v = exp(1j*2*pi*m*vd); 

        % Calculate target signal
        target_amplitude = sqrt(target_power);	% target amplitude (Swerling I model)
        target_signal = repmat(target_amplitude,num_samples,1).*v;	% target signal (Swerling I model)

        % Add the simulated target to the measured clutter
        H1 = H0 + target_signal;
        C_backup(t1:t2,(49 - RangeBin)) = H1;

    %% Insert second target (away from radar)
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

%% OSGO-CFAR on target dataset; for PD
kc1_mat_rvs = sort(kc1_mat,'descend');
kc2_mat_rvs = sort(kc2_mat,'descend');

% CFAR detection
for RangeBin = 1:NumOfRangeBins	% test a small sample to quicken processing
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
     
    %% Determine PD; consider all cases of targets crossing, not crossing, etc. (depends on speed)
        if (abs(target_doppler) <= 20)
            % check fD and fD1
            difference2 = signal_minus_T(fD1,kc1_mat(1):kc2_mat(1));
            difference1 = signal_minus_T(fD,(kc1_mat(48)):(kc2_mat(48)));
            target_detections = length(find(difference2 > 0));
            low_freq_1_1 = target_detections/length(kc1_mat(1):kc2_mat(1));
            target_detections = length(find(difference1 > 0));
            low_freq_48_1 = target_detections/length(kc1_mat(48):kc2_mat(48));
            % check fD+1 and fD1+1
            difference2 = signal_minus_T(fD1+1,kc1_mat(1):kc2_mat(1));
            difference1 = signal_minus_T(fD+1,(kc1_mat(48)):(kc2_mat(48)));
            target_detections = length(find(difference2 > 0));
            low_freq_1_2 = target_detections/length(kc1_mat(1):kc2_mat(1));
            target_detections = length(find(difference1 > 0));
            low_freq_48_2 = target_detections/length(kc1_mat(48):kc2_mat(48));
            % check fD+2 and fD1+2
            difference2 = signal_minus_T(fD1+2,kc1_mat(1):kc2_mat(1));
            difference1 = signal_minus_T(fD+2,(kc1_mat(48)):(kc2_mat(48)));
            target_detections = length(find(difference2 > 0));
            low_freq_1_3 = target_detections/length(kc1_mat(1):kc2_mat(1));
            target_detections = length(find(difference1 > 0));
            low_freq_48_3 = target_detections/length(kc1_mat(48):kc2_mat(48));
            % check fD-1 and fD1-1
            difference2 = signal_minus_T(fD1-1,kc1_mat(1):kc2_mat(1));
            difference1 = signal_minus_T(fD-1,(kc1_mat(48)):(kc2_mat(48)));
            target_detections = length(find(difference2 > 0));
            low_freq_1_4 = target_detections/length(kc1_mat(1):kc2_mat(1));
            target_detections = length(find(difference1 > 0));
            low_freq_48_4 = target_detections/length(kc1_mat(48):kc2_mat(48));
            % calculate best PD
            PD_mat1(1) = max([low_freq_1_1, low_freq_1_2, low_freq_1_3, low_freq_1_4]);
            PD_mat1(48) = max([low_freq_48_1, low_freq_48_2, low_freq_48_3, low_freq_48_4]);
        elseif fD == NFFT % target_doppler = 2500 Hz
            difference2 = signal_minus_T((fD1+1),kc1_mat(RangeBin):kc2_mat(RangeBin));
            difference1 = signal_minus_T(fD,(kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin)));
            target_detections = (length(find(difference1 > 0)) + length(find(difference2 > 0)));
            PD_mat1(RangeBin) = target_detections/(length(kc1_mat(RangeBin):kc2_mat(RangeBin)) + length((kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin))));
        % targets do not cross, do nothing where there is no data
        elseif (RangeBin > (RBMax)) && (RangeBin < (49 - (RBMax))) && (RBMax > 0) && (RBMax < 25)
            continue
        % targets do not cross, check RHS
        elseif (RangeBin <= (RBMax)) && (RBMax > 0) && (RBMax < 25)
            difference1 = signal_minus_T(fD1,kc1_mat(RangeBin):kc2_mat(RangeBin));
            difference3 = signal_minus_T((fD1+1),kc1_mat(RangeBin):kc2_mat(RangeBin));
            difference5 = signal_minus_T((fD1+2),kc1_mat(RangeBin):kc2_mat(RangeBin));
            difference7 = signal_minus_T((fD1-1),kc1_mat(RangeBin):kc2_mat(RangeBin));
            
            diff1_max = max([length(find(difference1 > 0)),length(find(difference3 > 0)),length(find(difference5 > 0)),length(find(difference7 > 0))]);
            target_detections = diff1_max;
            PD_mat1(RangeBin) = target_detections/(length(kc1_mat(RangeBin):kc2_mat(RangeBin)));
        % targets do not cross, check RHS
        elseif (RangeBin >= (49 - (RBMax))) && (RBMax > 0) && (RBMax < 25)
            difference2 = signal_minus_T(fD1,(kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin)));
            difference4 = signal_minus_T((fD1+1),(kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin)));
            difference6 = signal_minus_T((fD1+2),(kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin)));
            difference8 = signal_minus_T((fD1-1),(kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin)));
            
            diff2_max = max([length(find(difference2 > 0)),length(find(difference4 > 0)),length(find(difference6 > 0)),length(find(difference8 > 0))]);
            target_detections = diff2_max;
            PD_mat1(RangeBin) = target_detections/(length((kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin))));
        % targets cross, check RHS of no overlap
        elseif (RangeBin > (RBMax)) && (RBMax > 24)
            difference2 = signal_minus_T(fD,(kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin)));
            difference4 = signal_minus_T((fD+1),(kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin)));
            difference6 = signal_minus_T((fD+2),(kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin)));
            difference8 = signal_minus_T((fD-1),(kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin)));
            
            diff2_max = max([length(find(difference2 > 0)),length(find(difference4 > 0)),length(find(difference6 > 0)),length(find(difference8 > 0))]);
            target_detections = diff2_max;
            PD_mat1(RangeBin) = target_detections/(length((kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin))));
        % targets cross, check LHS of no overlap
        elseif (RangeBin < (49 - (RBMax))) && (RBMax > 24)
            difference1 = signal_minus_T(fD1,kc1_mat(RangeBin):kc2_mat(RangeBin));
            difference3 = signal_minus_T((fD1+1),kc1_mat(RangeBin):kc2_mat(RangeBin));
            difference5 = signal_minus_T((fD1+2),kc1_mat(RangeBin):kc2_mat(RangeBin));
            difference7 = signal_minus_T((fD1-1),kc1_mat(RangeBin):kc2_mat(RangeBin));
            
            diff1_max = max([length(find(difference1 > 0)),length(find(difference3 > 0)),length(find(difference5 > 0)),length(find(difference7 > 0))]);
            target_detections = diff1_max;
            PD_mat1(RangeBin) = target_detections/(length(kc1_mat(RangeBin):kc2_mat(RangeBin)));
        % targets cross, check overlaps
        else
            difference1 = signal_minus_T(fD1,kc1_mat(RangeBin):kc2_mat(RangeBin));
            difference2 = signal_minus_T(fD,(kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin)));
            difference3 = signal_minus_T((fD1+1),kc1_mat(RangeBin):kc2_mat(RangeBin));
            difference4 = signal_minus_T((fD+1),(kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin)));
            difference5 = signal_minus_T((fD1+2),kc1_mat(RangeBin):kc2_mat(RangeBin));
            difference6 = signal_minus_T((fD+2),(kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin)));
            difference7 = signal_minus_T((fD1-1),kc1_mat(RangeBin):kc2_mat(RangeBin));
            difference8 = signal_minus_T((fD-1),(kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin)));
            
            diff1_max = max([length(find(difference1 > 0)),length(find(difference3 > 0)),length(find(difference5 > 0)),length(find(difference7 > 0))]);
            diff2_max = max([length(find(difference2 > 0)),length(find(difference4 > 0)),length(find(difference6 > 0)),length(find(difference8 > 0))]);
            target_detections = diff1_max + diff2_max;
            PD_mat1(RangeBin) = target_detections/(length(kc1_mat(RangeBin):kc2_mat(RangeBin)) + length((kc1_mat(49-RangeBin)):(kc2_mat(49-RangeBin))));
        end
    
    clc;
    fprintf('Processing: %0.2f %% complete \n',RangeBin*100/NumOfRangeBins);
end

PD_mean = mean(PD_mat1(PD_mat1 > 0));   % average PD from all range bins

end