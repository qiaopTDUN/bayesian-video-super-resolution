addpath('mex');

% load the two frames
im1 = im2double(imread('fish1.jpg'));
im2 = im2double(imread('fish2.jpg'));
s   = 3;

im1lr = imresize(im1,1/s,'bicubic');
im2lr = imresize(im2,1/s,'bicubic');

% set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
alpha = 0.012;
ratio = 0.75;
minWidth = 20;
nOuterFPIterations = 7;
nInnerFPIterations = 1;
nSORIterations = 30;

para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];

% this is the core part of calling the mexed dll file for computing optical flow
% it also returns the time that is needed for two-frame estimation
tic;
[vx,vy,warpI2lr] = Coarse2FineTwoFrames(im1,im2,para);
toc

tic;
[vxlr,vylr,warpI2] = Coarse2FineTwoFrames(im1lr,im2lr,para);
toc

vx_hat = imresize(vxlr,size(vx));
vy_hat = imresize(vylr,size(vy));

dx     = vx - vx_hat;
dy     = vy - vy_hat;
d      = dx(:).^2+dy(:).^2;
mse    = sum(d(:))/numel(d)