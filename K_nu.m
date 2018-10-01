%% File information
% Written by F. Gini and M. S. Greco
% Modified by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% K_nu.m: Calculate the shape parameter nu of a K distribution by method of moments

function y = K_nu(x,m1,m2)
y = (4.*x.*gamma(x).^2)./(pi.*gamma(x + 0.5).^2) - m2./(m1.^2);
