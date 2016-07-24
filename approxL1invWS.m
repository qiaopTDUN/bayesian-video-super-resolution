function [W] = approxL1invWS(I, epsilon, hsz)
    dx  = [0 0 0; -1 1 0; 0 0 0]; % for x-derivative filter
    dy  = [0 -1 0; 0 1 0; 0 0 0]; % for y-derivative filter
    I   = reshape(I, hsz);
    DxI = conv2d(I, dx);
    DyI = conv2d(I, dy);
    md  = (DxI.^2 + DyI.^2 + epsilon^2) .^ (-0.5);
    np  = prod(size(md));
    W   = sparse(1:np, 1:np, md(:), np, np);
end