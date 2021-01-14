function x = tight_sspade(yc, masks, F, s, r, eps, max_iter)
%TIGHT_ASPADE Implentation of SSPADE declipping algorithm
%   yc : signal to declip
%   masks : structure describing reliable and clipped samples
%   s, r, eps : algorithm parameters
%   max_iter : maximum number of iterations
i = 1;
z0 = frana(F,yc);
zi = z0;
ui = zeros(size(zi));
k = s;
conv = false;
while i <= max_iter
    zb = hard_thresh(zi + ui, k);
   
    % projection
    v = frsyn(F,zb - ui);
    v(masks.Ir) = yc(masks.Ir);
    v(masks.Icm) = min(v(masks.Icm), yc(masks.Icm));
    v(masks.Icp) = max(v(masks.Icp), yc(masks.Icp));
    zi = zb - ui - frana(F, frsyn(F,zb - ui) - v);
    
    if norm(zi - zb) <= eps
        conv = true;
        break
    end
    ui = ui + zi - zb;
    i = i + 1;
    if mod(i,r) == 0
        k = k + s;
    end
end
x = frsyn(F,zi);

if ~conv
    disp('no convergence')
end
end

