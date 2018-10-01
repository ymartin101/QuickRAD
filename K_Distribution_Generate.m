%% File information
% Written by F. Gini and M. S. Greco
% Modified by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% K_Distribution_Generate.m: Generate clutter data from K distribution parameters, plot the PDF and estimate its generated parameters

%% Generate values from a K distribution
Nt = 1e6;
zc = (randn(Nt,1) + 1i*randn(Nt,1));
nu_K = 0.2;   % nu is the shape parameter
mu_K = 2.5;   % mu is the scale parameter
N = 16;       % time window

for n = 1:N:(Nt - N + 1)
    tau = gamrnd(nu_K,mu_K/nu_K,1,1); % texture is from a Gamma distribution
    z(n:(n + N - 1)) = sqrt(tau)*(zc(n:(n + N - 1))); 
end

%% PDF histogram
N = length(z);	% number of samples
y = abs(z);     % data is already K distributed
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
disp(['K distribution: nu = ' num2str(roundn(nu_K,-3)) ', nu_est = ' num2str(roundn(nu,-3)) ]);
disp(['K distribution: mu = ' num2str(roundn(mu_K,-3)) ', mu_est = ' num2str(roundn(mu,-3)) ]);
