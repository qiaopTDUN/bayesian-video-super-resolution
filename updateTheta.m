function [th] = updateTheta(I0, VX, VY, h2d, TP, params)
hsz   = TP.HSZ;
S     = TP.S;

Nq    = prod(TP.LSZ);
alpha = params.alpha;
beta  = params.beta;
th    = zeros(params.nfore + params.nback + 1, 1);
i     = -params.nfore:params.nback;

for j = 1:length(i)
    vx    = VX(:, j);
    vy    = VY(:, j);
    J     = TP.NSY(:, j);
    
    Fw    = generate_warp_matrix(vx, vy, hsz);
    FI    = reshape(Fw * I0(:), hsz);
    KFI   = conv2d(FI, h2d);
    SKFI  = S * KFI(:);
    d     = abs(J - SKFI);
    x_hat = sum(d(:));
    
    th(j) = (alpha + Nq - 1) / (beta + x_hat);
end
end