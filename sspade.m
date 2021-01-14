function x = sspade(yc, masks, F, s, r, eps, max_iter)
%SSPADE Summary of this function goes here
%   Detailed explanation goes here
N = length(yc);
D = frsynmatrix(F, N);

Mr = eye(N, N);
Mr = Mr(masks.Ir,:);
Mcp = eye(N, N);
Mcp = Mcp(masks.Icp,:);
Mcm = eye(N, N);
Mcm = Mcm(masks.Icm,:);

A = [-Mcp;Mcm]*D;
Ap = [real(A) -imag(A)];
b = [-Mcp;Mcm]*yc;
bp = b;
Aeq = Mr*D;
beq = Mr*yc;
Aeqp = [real(Aeq) -imag(Aeq); imag(Aeq) real(Aeq)];
beqp = [beq; zeros(size(beq))];

i = 1;
z0 = frana(F,yc);
zi = z0;
ui = zeros(size(zi));

k = s;

conv = false;
while i <= max_iter
    zb = hard_thresh(zi + ui, k);
   
    [zi,~,~,exitflag] = lsqlin(eye(2*length(zb)),[real(zb); imag(zb)] - [real(ui); imag(ui)],Ap,bp,Aeqp,beqp,[],[],[real(zi); imag(zi)],optimoptions(@lsqlin,'Algorithm', 'interior-point', 'OptimalityTolerance', 1e-4, 'MaxIterations', 10, 'Display','off'));
    if exitflag ~= 1
        disp(['exit flag ' num2str(exitflag)])
    end
    zi = zi(1:length(zi)/2) + 1j*zi(length(zi)/2+1:end);
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
disp(conv)
end

