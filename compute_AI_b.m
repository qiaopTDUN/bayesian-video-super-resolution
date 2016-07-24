function [AI, b] = compute_AI_b(I0, VX, VY, TH, h2d, TP, params, Wi, Ws)
hsz  = TP.HSZ;
dx   = [0 0 0; -1 1 0; 0 0 0]; % for x-derivative filter
dy   = [0 -1 0; 0 1 0; 0 0 0]; % for y-derivative filter
S    = TP.S;

AI   = 0;
b    = 0;
i    = -params.nfore:params.nback;

for j = 1:length(i)
    vx = VX(:, j);
    vy = VY(:, j);
    J  = TP.NSY(:, j);
    th = TH(j);
    %% AI
    wi = Wi{j};
    Fw = generate_warp_matrix(vx, vy, hsz);
    FI = reshape(Fw * I0(:), hsz);
    KFI = conv2d(FI, h2d);
    SKFI = S * KFI(:);
    WiSKFI = wi * SKFI;
    StWiSKFI = reshape(S' * WiSKFI, hsz);
    KtStWiSKFI = conv2dt(StWiSKFI, h2d);
    thFwtKtStWiSKFI = reshape(th * Fw' * KtStWiSKFI(:), hsz);
    %% b
    WiJi = wi * J(:);
    StWiJi = reshape(S' * WiJi, hsz);
    KtStWiJi = conv2dt(StWiJi, h2d);
    thFwtKtStWiJi = reshape(th * Fw' * KtStWiJi(:), hsz);
    %% accu
    AI = AI + thFwtKtStWiSKFI;
    b  = b + thFwtKtStWiJi;
end

DxI = conv2d(reshape(I0, hsz), dx);
WsDxI = reshape(Ws * DxI(:), hsz);
DxtWsDxI = conv2dt(WsDxI, dx);
DyI = conv2d(reshape(I0, hsz), dy);
WsDyI = reshape(Ws * DyI(:), hsz);
DytWsDyI = conv2dt(WsDyI, dy);
AI = AI + params.ita * (DxtWsDxI + DytWsDyI);
end