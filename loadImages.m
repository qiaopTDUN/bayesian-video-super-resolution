function [TP] = loadImages(params)
gt    = [];
input = [];
noisy = [];
ofx   = [];
ofy   = [];

hsz   = [];
lsz   = [];
nsz   = [];

scl   = params.scaleup;

h_2D_gt = fspecial('gaussian', [params.hsize params.hsize], params.hsigma);
h_2D    = fspecial('gaussian', [params.hsize params.hsize], params.hsigma_in);
%% parameters for celiu optical flow
alpha = 0.012;
ratio = 0.75;
minWidth = 20;
nOuterFPIterations = 7;
nInnerFPIterations = 1;
nSORIterations = 30;
para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];

for i = -params.temporalR:params.temporalR
    path   = sprintf('%s %03d.png', params.fn, i + params.idx0);
    I0     = imread(path);
    if size(I0,3) > 1
        I1 = rgb2ycbcr(I0);
        HR = im2double(I1(:,:,1)); 
    else
        HR = im2double(I0);
    end
    
    hsz    = size(HR);
    gt     = [gt HR(:)];
end

S = generate_downsample_matrix(scl, hsz);

for i = 1:2 * params.temporalR + 1
    HR     = reshape(gt(:,i), hsz);
    KI     = conv2d(HR, h_2D_gt);
    LR     = reshape(S * KI(:), hsz/scl);
    INTERP = imresize(LR, scl, 'bicubic');

    lsz    = size(LR);
    nsz    = size(INTERP);
    
    noisy  = [noisy LR(:)];
    input  = [input INTERP(:)];
end

referI = reshape(noisy(:, params.idx0), lsz);
theta  = zeros(params.temporalR*2+1, 1);
Nq     = prod(hsz);

for i = -params.temporalR:params.temporalR
    if i==0
        tx  = zeros(hsz);
        ty  = tx;
        ofx = [ofx tx(:)];
        ofy = [ofy ty(:)];
        nextI = reshape(noisy(:, params.idx0 + i + 1), lsz);
        lastI = reshape(noisy(:, params.idx0 + i - 1), lsz);
        d     = (lastI - nextI) / 2;
        theta(i+params.idx0) = Nq / sum(abs(d(:)));
        continue;
    end
    warpI = reshape(noisy(:, params.idx0 + i), lsz);
    [tx, ty, ~] = Coarse2FineTwoFrames(warpI, referI, para);
    
    tx  = imresize(tx, params.scaleup, 'bicubic');
    ty  = imresize(ty, params.scaleup, 'bicubic');
    ofx = [ofx tx(:)];
    ofy = [ofy ty(:)];
    
    if i<0
        nextI = reshape(noisy(:, params.idx0 + i + 1), lsz);
        d     = warpI - nextI;
        theta(i+params.idx0) = Nq / sum(abs(d(:)));
    elseif i>0
        lastI = reshape(noisy(:, params.idx0 + i - 1), lsz);
        d     = lastI - warpI;
        theta(i+params.idx0) = Nq / sum(abs(d(:)));
    end
end

TP.GND_TRU   = gt;
TP.INPUT     = input;
TP.NOISY     = noisy;
TP.OPT_FLO_X = ofx;
TP.OPT_FLO_Y = ofy;
TP.THETA     = abs(theta);
TP.h_2D      = h_2D;
TP.h_2D_gt   = h_2D_gt;

TP.HSZ       = hsz;
TP.LSZ       = lsz;
TP.NSZ       = nsz;
TP.S         = S;
%% input
TP.I0        = input(:,params.idx0);
%% output
TP.Ik        = TP.I0;
TP.NSY       = noisy(:,params.idx0-params.nfore:params.idx0+params.nback);
TP.OFX       = ofx(:,params.idx0-params.nfore:params.idx0+params.nback);
TP.OFY       = ofy(:,params.idx0-params.nfore:params.idx0+params.nback);
TP.TH        = theta(params.idx0-params.nfore:params.idx0+params.nback);
TP.H2D       = h_2D;

TP.It        = cell(params.maxit_outer*2,1);
TP.Ht        = cell(params.maxit_outer*2,1);
end