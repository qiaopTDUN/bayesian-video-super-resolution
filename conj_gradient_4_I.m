function [Ik] = conj_gradient_4_I(AI, I0, b, Wi, Ws, vx, vy, th, h2d, TP, params)
%A = th(j)K'S'W0SK + eta*(Dx'WsDx + Dy'WsDy) + sum(th(j)F'K'S'W0SKF)
%b = th(j)K'S'W0J0;
%I = A \ b;

maxit = params.maxit_conjI;
eps   = params.epsilon_conjI;

r     = AI(:) - b(:);
p     = -r;
k     = 0;
I     = I0(:);
rsize = length(r);

while k < maxit
    [Ap,~]= compute_AI_b(p, vx, vy, th, h2d, TP, params, Wi, Ws);
    alpha = (r'*r) / (p'*Ap(:));
    I     = I + alpha * p;
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

Ik = I;
end