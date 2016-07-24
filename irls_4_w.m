function [ux, uy] = irls_4_w(I0, VX, VY, TH, h2d, TP, params)
hsz     = TP.HSZ;
np      = prod(hsz);

dx      = [0 0 0; -1 1 0; 0 0 0]; % for x-derivative filter
dy      = [0 -1 0; 0 1 0; 0 0 0]; % for y-derivative filter

% Ix      = conv2d(reshape(I0,hsz), dx);
% Iy      = conv2d(reshape(I0,hsz), dy);

maxit   = params.maxit_innerW;
eps     = params.epsilon_innerW;

ux      = VX;
uy      = VY;
noisy   = TP.NSY;
lambda  = params.lambda;

i = -params.nfore:params.nback;
for j = 1:length(i)
    if i(j) == 0
        continue;
    end
    
    vx = VX(:,j);
    vy = VY(:,j);
    J  = noisy(:,j);
    th = TH(j);
    zeta = lambda / th;
    
    Fwi = generate_warp_matrix(vx, vy, hsz);
    Ix  = conv2d(reshape(I0,hsz), dx);
    Iy  = conv2d(reshape(I0,hsz), dy);
    FIx = sparse(1:np, 1:np, Fwi * Ix(:), np, np);
    FIy = sparse(1:np, 1:np, Fwi * Iy(:), np, np);
    
    d0  = zeros(2*np,1);
    m   = 1;
    
    while m <= maxit
        dk_oldm = d0;
        [Wi, Wus, Wvs] = compute_Wi_Wus_Wvs(I0, J, vx, vy, h2d, TP, params);
        [Aw, b]        = compute_Aw_b(I0, J, d0, Fwi, FIx, FIy, vx, vy, zeta, h2d, TP, Wi, Wus, Wvs);
        d0             = conj_gradient_4_w(Aw, b, I0, J, d0, Fwi, FIx, FIy, vx, vy, zeta, h2d, Wi, Wus, Wvs, TP, params);
        
        vx = vx + d0(1:np);
        vy = vy + d0(np+1:2*np);
        
        rms = norm(d0(:) - dk_oldm(:)) / norm(dk_oldm(:));
        fprintf('%03d: \t%.6f\n', m, rms);

        if rms <= eps
            break;
        end

%         pause(1);
        m = m + 1;
    end
    ux(:, j) = vx;
    uy(:, j) = vy;
end
end