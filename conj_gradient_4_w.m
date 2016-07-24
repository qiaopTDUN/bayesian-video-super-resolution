function [dk] = conj_gradient_4_w(Aw, b, I0, J, d0, Fwi, FIx, FIy, vx, vy, zeta, h2d, Wi, Wus, Wvs, TP, params)
%    / Ix'K'S'WiSKIx + zeta(i)Lu   Ix'K'S'WiSKIy \
%A = |                                           |
%    \ Iy'K'S'WiSKIx   Iy'K'S'WiSKIy + zeta(i)Lv /
%b = th(j)K'S'W0J0;
%I = A \ b;

maxit = params.maxit_conjW;
eps   = params.epsilon_conjW;

r     = Aw(:) - b(:);
p     = -r;
k     = 0;
d     = d0(:);
rsize = length(r);

while k < maxit
    [Ap,~]= compute_Aw_b(I0, J, p, Fwi, FIx, FIy, vx, vy, zeta, h2d, TP, Wi, Wus, Wvs);
    alpha = (r'*r) / (p'*Ap(:));
    d     = d + alpha * p;
    r_new = r + alpha * Ap(:);
    beta  = (r_new' * r_new) / (r' * r);
    p_new = -r_new + beta * p;
    
    diff  = norm(r_new) / rsize;
    
    if diff < eps
        break;
    end
    
    k = k + 1;
    r = r_new;
    p = p_new;
end

dk = d;
end