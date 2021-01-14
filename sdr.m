function s = sdr(x, y, masks)
%SDR Summary of this function goes here
%   Detailed explanation goes here
Ic = masks.Icp | masks.Icm;
s = 20*log10(norm(x(Ic)) / norm(x(Ic) - y(Ic)));
end

