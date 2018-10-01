%% File information
% Written by F. Gini and M. S. Greco
% Modified by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% K_Distribution_Fit.m: Fit measured clutter data to a K distribution parameters, plot the PDF and return its parameters

%% Setup
load 'CFA17_002.mat';       % load dataset with variable 'Cdata' (clutter vector)
RangeBin = 1;               % select range bin to fit K distribution to

%% PDF histogram
N = length(Cdata(:,RangeBin));     % number of samples
y = abs(Cdata(:,RangeBin)).^2;     % square law
delta = 0.05;
yapl = 0:delta:max(y);
histo = hist(y,yapl)./(N*delta);

%% Moments
m(1) = mean(y);
mnorm = m(1);
m(2) = mean(y.^2);
i = 1:6;

%% K fitting
nu = fzero('K_nu',1,[],m(1),m(2));      % shape parameter of the clutter
mu = m(2)./2;                           % scale parameter of the clutter
yk = K_pdf(yapl,nu,mu);                 % returns K-PDF of the clutter
% K theoretical moments
mk = (2.*mu).^(i./2).*gamma(nu + i./2).*gamma(i./2 + 1)./(nu.^(i./2).*gamma(nu));

%% Plot PDF fitting
figure(1);
semilogy(yapl,histo,'*b',yapl,yk,'-.r','LineWidth',2);
axis([0 6 0.001 10]);
xlabel('Return Strength');
ylabel('Logarithmic Probability Density');
title('PDF Histogram and K-Distribution Fit');
legend('PDF Histogram','K-Distribution');
grid on;

%% Outputs
clc;
disp(['K distribution: nu = ' num2str(roundn(nu,-3)) ]);
disp(['K distribution: mu = ' num2str(roundn(mu,-3)) ]);
