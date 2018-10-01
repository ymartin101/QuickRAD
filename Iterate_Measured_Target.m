%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% Iterate_Measured_Target.m: Determine PFA for OSGO-CFAR for a measured target dataset for a varying parameter

%% PD plots for OS- and OSGO-CFAR for a "multiple" dataset
PFA_array = [];
variable_array = 8:8:32;        % change according to Variable under test, eg. N

for Variable = variable_array   % change Variable in "OSGO_Iterate" file
    PFA = OSGO_Iterate(Variable);
    PFA_array = [PFA_array PFA];
end
