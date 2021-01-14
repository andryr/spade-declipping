function [v] = hard_thresh(u,k)
%HARD_THRESH Summary of this function goes here
%   Detailed explanation goes here
l = flip(sort(abs(u)));
alpha = l(k);
v = u .* (abs(u) >= alpha);
end

