function [TPk] = mfsr_irls(TP, params)
esp   = params.epsilon_outer;
maxit = params.maxit_outer;
k     = 1;
I0    = TP.Ik;
vx    = TP.OFX;
vy    = TP.OFY;
th    = TP.TH;
h2d   = TP.H2D;

while k <= maxit
    fprintf('Iter: %d\n', k);
    imwrite(reshape(I0, TP.HSZ), ['estimate-' num2str(k) '.png']);
    TP.It{k} = reshape(I0, TP.HSZ);
    TP.Ht{k} = h2d;
    
    Ik_oldk  = I0;
    
    tic
    %% estimate optical flow
    fprintf('motion estimation...\n');
    [vx, vy] = irls_4_w(I0, vx, vy, th, h2d, TP, params);
    %% estimate noise level
    fprintf('noise level estimation...\n');
    th       = updateTheta(I0, vx, vy, h2d, TP, params);
    %% estimate image
    fprintf('image estimation...\n');
    I0       = irls_4_I(I0, vx, vy, th, h2d, TP, params);
    %% estimate kernel
    fprintf('kernel estimation...\n')
    h2d      = irls_4_k_full(I0, th, h2d, TP, params);
    %% check convergence
    rms      = norm(I0(:) - Ik_oldk(:)) / norm(Ik_oldk(:));
    T        = toc;
    fprintf('RMS: %.6f\tRuntime: %.0f\n', rms, T);

    if rms <= esp
        break;
    end
    
    k = k + 1;
end

imwrite(reshape(I0, TP.HSZ), ['estimate-' num2str(k) '.png']);
TP.It{k} = reshape(I0, TP.HSZ);
TP.Ht{k} = h2d;
TP.Ik    = I0;
TP.OFX   = vx; 
TP.OFY   = vy; 
TP.TH    = th; 
TP.H2D   = h2d;

TPk      = TP;
end