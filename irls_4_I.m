function [Ik] = irls_4_I(I0, vx, vy, th, h2d, TP, params)
m       = 1;
maxit   = params.maxit_innerI;
eps     = params.epsilon_innerI;

while m <= maxit
    Ik_oldm  = I0;
    
    [Wi, Ws] = compute_Wi_Ws(I0, vx, vy, h2d, TP, params);
    [AI, b]  = compute_AI_b(I0, vx, vy, th, h2d, TP, params, Wi, Ws);
    I0       = conj_gradient_4_I(AI, I0, b, Wi, Ws, vx, vy, th, h2d, TP, params);

    rms = norm(I0(:) - Ik_oldm(:)) / norm(Ik_oldm(:));
    fprintf('%03d: \t%.6f\n', m, rms);

    if rms <= eps
        break;
    end

    pause(1);
    m = m + 1;
end  
Ik = I0;
end