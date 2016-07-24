function [hk] = irls_4_k_full(I0, th, h0, TP, params)
m       = 1;
maxit   = params.maxit_innerK;
esp     = params.epsilon_innerK;    

hsz     = TP.HSZ;
h00     = otf2psf(psf2otf(h0, hsz));
while m <= maxit
    hk_oldm   = h00;
    [W0, Wks] = compute_W0_Wks_full(I0, h00, TP, params);
    [Ak, b]   = compute_Ak_b_full(I0, th, h00, TP, params, W0, Wks);
    h00       = conj_gradient_4_k_full(Ak, b, h00, I0, W0, Wks, th, TP, params);
    
    rms = norm(h00(:) - hk_oldm(:)) / norm(hk_oldm(:));
    fprintf('%03d: \t%.6f\n', m, rms);
%     h0(h0 < 0) = 0;
%     h0 = h0 / sum(h0);
    
    if rms <= esp
        break;
    end

    pause(1);
    m = m + 1;
end
half = (params.hsize-1)/2;
h00  = reshape(h00, hsz);
hk   = h00(hsz(1)/2-half + 1:hsz(1)/2+half + 1,hsz(2)/2-half + 1:hsz(2)/2+half + 1);
hk   = hk / sum(hk(:));
end