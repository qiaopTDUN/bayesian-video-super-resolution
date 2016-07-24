function [W0, Wks] = compute_W0_Wks_full(I0, h0, TP, params)
    KI  = conv2d(reshape(I0, TP.HSZ), reshape(h0, TP.HSZ));
    SKI = TP.S * KI(:);
    d   = SKI - TP.NOISY(:, params.idx0);
    md  = (d.^2 + params.epsilon^2) .^ (-0.5);
    np  = prod(size(md));
    W0  = sparse(1:np, 1:np, md, np, np);
    
    Wks = approxL1invWS(h0, params.epsilon, TP.HSZ); 
end