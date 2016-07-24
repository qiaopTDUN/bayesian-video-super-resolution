function [Aw, b] = compute_Aw_b(I0, J, d0, Fwi, FIx, FIy, vx, vy, zeta, h2d, TP, Wi, Wus, Wvs)
hsz  = TP.HSZ;
np   = prod(hsz);
dx   = [0 0 0; -1 1 0; 0 0 0]; % for x-derivative filter
dy   = [0 -1 0; 0 1 0; 0 0 0]; % for y-derivative filter
S    = TP.S;

du   = d0(1:np);
dv   = d0(np+1:2*np);
%% Aw
A111 = conv2d(reshape(FIx * du, hsz), h2d);   % K * FIx * Du
A111 = conv2dt(reshape(S' * Wi * S * A111(:), hsz), h2d);    % K' * S' * Wi * S * KFIxDu
A111 = FIx' * A111(:);    % FIx' * K'S'WiSKFIxDu
A112 = conv2d(reshape(du, hsz), dx);    % Dx * Du
A112 = zeta * conv2dt(reshape(Wus * A112(:), hsz), dx);    % zeta * Dx' * Wus * DxDu;
A113 = conv2d(reshape(du, hsz), dy);    % Dy * Du
A113 = zeta * conv2dt(reshape(Wus * A113(:), hsz), dy);    % zeta * Dy' * Wus * DyDu;
A11  = A111(:) + A112(:) + A113(:);

A12  = conv2d(reshape(FIy * dv, hsz), h2d);    % K * FIy * Dv
A12  = conv2dt(reshape(S' * Wi * S * A12(:), hsz), h2d);    % K' * S' Wi * S * KFIyDv
A12  = FIx' * A12(:);    % FIx' * K'S'WiSKFIyDv

A21  = conv2d(reshape(FIx * du, hsz), h2d);    % K * FIx * Du
A21  = conv2dt(reshape(S' * Wi * S * A21(:), hsz), h2d);    % K' * S' Wi * S * KFIxDu
A21  = FIy' * A21(:);    % FIy' * K'S'WiSKFIxDu

A221 = conv2d(reshape(FIy * dv, hsz), h2d);   % K * FIy * Dv
A221 = conv2dt(reshape(S' * Wi * S * A221(:), hsz), h2d);    % K' * S' * Wi * S * KFIyDv
A221 = FIy' * A221(:);    % FIy' * K'S'WiSKFIyDv
A222 = conv2d(reshape(dv, hsz), dx);    % Dx * Dv
A222 = zeta * conv2dt(reshape(Wvs * A222(:), hsz), dx);    % zeta * Dx' * Wvs * DxDv;
A223 = conv2d(reshape(dv, hsz), dy);    % Dy * Dv
A223 = zeta * conv2dt(reshape(Wvs * A223(:), hsz), dy);    % zeta * Dy' * Wvs * DyDv;
A22  = A221(:) + A222(:) + A223(:);

A1   = A11 + A12;
A2   = A21 + A22;
%% b
b111 = conv2d(reshape(vx, hsz), dx);    % Dx * vx
b111 = zeta * conv2dt(reshape(Wus * b111(:), hsz), dx);    % zeta * Dx' * Wus * Dxvx
b112 = conv2d(reshape(vx, hsz), dy);    % Dy * vx
b112 = zeta * conv2dt(reshape(Wus * b112(:), hsz), dy);    % zeta * Dy' * Wus * Dyvx
b11  = b111(:) + b112(:);

b211 = conv2d(reshape(vy, hsz), dx);    % Dx * vy
b211 = zeta * conv2dt(reshape(Wvs * b211(:), hsz), dx);    % zeta * Dx' * Wus * Dxvy
b212 = conv2d(reshape(vy, hsz), dy);    % Dy * vy
b212 = zeta * conv2dt(reshape(Wvs * b212(:), hsz), dy);    % zeta * Dy' * Wus * Dyvy
b21  = b211(:) + b212(:);

c1   = conv2d(reshape(Fwi * I0, hsz), h2d);    % K * Fwi * I
c1   = S' * Wi * S * c1(:);    % S' * Wi * S * KFwiI
c1   = conv2dt(reshape(c1, hsz), h2d);    % K' * S'WiSKFwiI
c2   = S' * Wi * J(:);    % S' * Wi * J
c2   = conv2dt(reshape(c2, hsz), h2d);    % K' * S'WiJ
c    = c1(:) - c2(:);
b12  = FIx' * c;
b22  = FIy' * c;
b1   = - b11 - b12;
b2   = - b21 - b22;

%%
Aw   = [A1;A2];
b    = [b1;b2];
end