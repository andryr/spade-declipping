function x = tight_aspade(yc, masks, F, s, r, eps, max_iter)
%TIGHT_ASPADE Implentation of ASPADE declipping algorithm
%   yc : signal to declip
%   masks : structure describing reliable and clipped samples
%   s, r, eps : algorithm parameters
%   max_iter : maximum number of iterations
i = 1;
xi = yc;
Axi = frana(F, xi);
ui = zeros(size(Axi));
k = s;
conv = false;
while i <= max_iter
    zi = hard_thresh(Axi + ui, k);

    % projection
    v = frsyn(F,zi - ui);
    v(masks.Ir) = yc(masks.Ir);
    v(masks.Icm) = min(v(masks.Icm), yc(masks.Icm));
    v(masks.Icp) = max(v(masks.Icp), yc(masks.Icp));
    xi = v;
    
    Axi = frana(F,xi);
    if norm(Axi - zi) <= eps
        %disp(['converged in ' num2str(i) ' iterations'])
        conv = true;
        break
    end

    ui = ui + Axi - zi;
    i = i + 1;
    if mod(i,r) == 0
        k = k + s;
    end
end
x = xi;

if ~conv
    disp('no convergence')
end
end

