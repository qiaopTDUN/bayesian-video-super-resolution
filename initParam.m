function [params] = initParam()
    %  model parameters
    params.fn        = 'data/calendar/original/Frame';
    params.scaleup   = 2;
    params.idx0      = 16;
    params.temporalR = 15;
    params.nfore     = 2;  %15
    params.nback     = 2;  %15
    params.ita       = 0.02;  %0.02
    params.lambda    = 1;  %1
    params.ksi       = 0.7;  %0.7
    params.alpha     = 1;  %1
    params.beta      = 0.1;  %0.1
    params.hsize     = 15;
    params.hsigma    = 1.6;
    params.hsigma_in = 2;
    params.epsilon   = 0.001;
    % optimization stop criterion
    params.epsilon_outer  = 1e-4;
    params.epsilon_innerI = 5e-4;
    params.epsilon_innerW = 1e-4;
    params.epsilon_innerK = 1e-4;
    params.epsilon_conjI  = 1e-1;
    params.epsilon_conjW  = 1e-3;
    params.epsilon_conjK  = 1e-2;
    params.maxit_outer    = 50;
    params.maxit_innerI   = 5;
    params.maxit_innerW   = 3;
    params.maxit_innerK   = 3;
    params.maxit_conjI    = 5;
    params.maxit_conjW    = 5;
    params.maxit_conjK    = 100;
end