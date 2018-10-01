%% File information
% Written by F. Gini and M. S. Greco
% Modified by M. Y. Martin (MRTMOG014)
% EEE4022S (2018)
% K_pdf.m: Return the K pdf values using parameters nu and mu at the values in x

function y = K_pdf(x,nu,mu)
%   Returns the K pdf with parameters nu and mu at the values in x
%   The size of P is the common size of the input arguments
%   A scalar input functions as a constant matrix of the same size as the other inputs

if nargin < 3
    error('Requires three input arguments.'); 
end

[errorcode,x,nu,mu] = distchck(3,x,nu,mu);
if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

y = zeros(size(x));

k1 = find(nu <= 0 | mu <= 0);
if any(k1)
   tmp   = NaN;
   y(k1) = tmp(ones(size(k1)));
end

k = find(x > 0 & nu > 0 & mu > 0);
if any(k)
   y(k) = sqrt(2.*nu(k)./mu(k))./(gamma(nu(k)).*2.^(nu(k) - 1)).*...
  (sqrt(2.*nu(k)./mu(k)).*x(k)).^nu(k).*besselk(nu(k) - 1,(sqrt(2.*nu(k)./mu(k)).*x(k)));
end
