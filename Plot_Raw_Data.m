%% File information
% Written by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% Plot_Raw_Data.m: Produce a time and range bin plot of a measured dataset

%% Startup
close all;
clc;
% select dataset to plot
load('TFC15_008.mat');  % loads variables: Cdata, NumOfPRIs, NumOfRangeBins, PRI_s; must be included in .mat file

%% Plot raw measured data
range_bins = 1:1:size(Cdata,2);
time = (1:1:size(Cdata,1))*PRI;
figure;
imagesc(range_bins,time,20*log10(abs(Cdata)));
colorbar;
colormap(jet);
hold on;
title('Measured Data');
xlabel('Range Bin');
ylabel('Time [seconds]');
