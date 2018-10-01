%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% PD_vs_Speed_Multiple.m: Plot PD vs Target Speed for OS- and OSGO-CFAR for a "multiple" dataset

%% PD plots for OS- and OSGO-CFAR for a "multiple" dataset
PD_array_OS = [];
PD_array_OSGO = [];
variable_array = 0:50:2500;	% change according to Variable under test

for Variable = variable_array   % change Variable to target_doppler in both function files
    PD = PD_Multiple_OS(Variable);
    PD_array_OS = [PD_array_OS PD];
    PD = PD_Multiple(Variable);
    PD_array_OSGO = [PD_array_OSGO PD];
end

c0 = 299792458;                             % speed of light
lambda = c0/(6.9*(10^9));                   % wavelength of the radar; f = 6.9 GHz

%% Plot PD vs Target Speed for OS- and OSGO-CFAR
figure;
variable_array = variable_array.*lambda./2;
plot(variable_array,PD_array_OSGO,variable_array,PD_array_OS,'LineWidth',1.5);
title('Probability of Detection vs Target Speed');
xlabel('Target Speed [m/s]');
PD_vs_speed_legend = legend('OSGO-CFAR','OS-CFAR');
set(PD_vs_speed_legend,'Location','SouthEast');
ylabel('PD');
grid on;
ylim([0 1]);
