function y = redundant_dft(x, r)
%REDUNDANT_DFT Summary of this function goes here
%   Detailed explanation goes here
n = length(x);
nr = n*r;
z = zeros([nr 1]);
z(1:n) = x;
y = fft(z) / sqrt(nr);
end

