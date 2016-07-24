function [Ak, b] = compute_Ak_b_full(I0, th, h0, TP, params, W0, Wks)
dx  = [0 0 0; -1 1 0; 0 0 0]; % for x-derivative filter
dy  = [0 -1 0; 0 1 0; 0 0 0]; % for y-derivative filter

hsz = TP.HSZ;
S   = TP.S;
j   = params.idx0;
i   = -params.nfore:params.nback;
th0 = th(i==0);
J   = TP.NOISY(:, j);

%% Ak
KI  = conv2d(reshape(I0, hsz), reshape(h0, hsz));
StW0SKI = S' * W0 * S * KI(:);
UtStW0SKI = th0 * conv2dt(reshape(StW0SKI, hsz), reshape(I0, hsz));

Dxk  = conv2d(reshape(h0, hsz), dx);
WksDxk = Wks * Dxk(:);
DxtWksDxk = params.ksi * conv2dt(reshape(WksDxk, hsz), dx);
Dyk  = conv2d(reshape(h0, hsz), dy);
WksDyk = Wks * Dyk(:);
DytWksDyk = params.ksi * conv2dt(reshape(WksDyk, hsz), dy);
Ak = UtStW0SKI(:) + DxtWksDxk(:) + DytWksDyk(:);
Ak = UtStW0SKI(:) + params.ksi * ones(size(UtStW0SKI(:)));
%% b
StW0J = S' * W0 * J;
b = th0 * conv2dt(reshape(StW0J, hsz), reshape(I0, hsz));
end