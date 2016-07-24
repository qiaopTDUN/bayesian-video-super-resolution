function [W] = approxL1inv(I, J, tx, ty, S, h2d, epsilon, hsz)
    Fwi  = generate_warp_matrix(tx, ty, hsz);
    FI   = reshape(Fwi * I(:), hsz);
    KFI  = conv2d(FI, h2d);
    SKFI = S * KFI(:);
    d    = SKFI - J(:);
    md   = (d.^2 + epsilon^2) .^ (-0.5);
    np   = prod(size(md));
    W    = sparse(1:np, 1:np, md, np, np);
end