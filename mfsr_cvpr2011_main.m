clear

addpath(genpath('C:\Codes\mfsr-cvpr2011-repo\celiu_optical_flow\'));
%% parameters initialization
[params] = initParam();

%% load data and initialize
[TP] = loadImages(params);

%% optimize the model parameters via IRLS
[TP] = mfsr_irls(TP, params);

save 'model.mat' 'TP' 'params'