function x = redundant_idft(y, r)
%REDUNDANT_IDFT Summary of this function goes here
%   Detailed explanation goes here
nr = length(y);
n = nr / r;
x = ifft(y);
x = x(1:n) * sqrt(nr);
end

