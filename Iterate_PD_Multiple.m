%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% Iterate_PD_Multiple.m: Iterate through PD_Multiple to find PD and PFA for a varying parameter

%% Iterate through PD_Multiple
PD_array = [];
PFA_array = [];
variable_array = 8:2:32;        % change according to Variable under test, eg. N

for Variable = variable_array   % change Variable in "PD_Multiple"" file
    [PD,PFA] = PD_Multiple(Variable);
    % add to PD and PFA array to track them as the parameter varies
    PD_array = [PD_array PD];
    PFA_array = [PFA_array PFA];
end
