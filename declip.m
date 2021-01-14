function [x, masks] = declip(yc, params)
%DECLIP Summary of this function goes here
%   Detailed explanation goes here
tau = max(abs(yc));
masks.Icp = yc == tau;
masks.Icm = yc == -tau;
masks.Ir = ~(masks.Icp | masks.Icm);

M = params.window_length;
a = params.hop_length;
L = length(yc);
N = L/a;
x = zeros(size(yc));

win = sqrt(hamming(M));

if strcmp(params.algo, 'aspade')
    spade = @tight_aspade;
elseif strcmp(params.algo, 'sspade')
    spade = @tight_sspade;
end
norm_win = zeros(size(x));
for i=0:N-1
    idx1 = i*a + 1;
    idx2 = min(i*a + M, L);
    yc_chunk = zeros([M 1]);
    yc_chunk(1:min(M,L-i*a)) = yc(idx1:idx2);
    yc_chunk = yc_chunk.*win;
    
    chunk_masks.Ir = masks.Ir(idx1:idx2);
    chunk_masks.Icp = masks.Icp(idx1:idx2);
    chunk_masks.Icm = masks.Icm(idx1:idx2);
    x_chunk = spade(yc_chunk, chunk_masks, params.frame, params.s, params.r, params.eps, params.max_iter);
    x(idx1:idx2) = x(idx1:idx2) + x_chunk(1:min(M,L-i*a)) .* win((1:min(M,L-i*a)));
    norm_win(idx1:idx2) = norm_win(idx1:idx2) + win(1:min(M,L-i*a)).^2;
end
x = x ./ norm_win;
end

