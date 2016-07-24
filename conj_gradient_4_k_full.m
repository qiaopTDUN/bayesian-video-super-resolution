function [hk] = conj_gradient_4_k_full(Ak, b, h0, I0, W0, Wks, th, TP, params)
%A = th(j)U'S'W0SU + ksi*(Dx'WsDx + Dy'WsDy)
%b = th(j)U'S'W0J0;
%I = A \ b;

maxit = params.maxit_conjK;
eps   = params.epsilon_conjK;

r     = Ak(:) - b(:);
p     = -r;
k     = 0;
h     = h0(:);
rsize = length(r);

while k < maxit
    [Ap,~]= compute_Ak_b_full(I0, th, p, TP, params, W0, Wks);
    alpha = (r'*r) / (p'*Ap(:));
    h     = h + alpha * p;
    r_new = r + alpha * Ap(:);
    beta  = (r_new' * r_new) / (r' * r);
    p_new = -r_new + beta * p;
    
    h(h < 0) = 0;
    h     = h / sum(h);
    
    diff  = norm(r_new) / rsize;
    
    if diff < eps
        break;
    end
    
    k = k + 1;
    r = r_new;
    p = p_new;
end

hk = h;
end