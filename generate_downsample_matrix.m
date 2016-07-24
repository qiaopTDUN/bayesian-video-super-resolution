function [S] = generate_downsample_matrix(scl, hsz)
    slect  = zeros(hsz);
    slect(1:scl:end,1:scl:end,:) = 1;
    imdims = size(slect);
    I      = speye(prod(imdims),prod(imdims));
    sub    = find(slect(:) == 1);
    S      = I(sub,:);
end